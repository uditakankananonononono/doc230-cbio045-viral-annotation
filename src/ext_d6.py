import json
SRC=open('src/run.py').read()
exec(SRC.split("net=Net(); opt")[0])   # setup identical to run.py: data, splits, B1, windows(), kfeat, lr, Net
mid=SRC.split("def predict_record")[1].split("for name,rs in (('d2_cold'")[0]
exec("def predict_record"+mid)
def metr(rs):
    Y=[];M=[]
    for r in rs:
        o=predict_record(r)
        if o is None: continue
        a,b,ok=o; Y.append(r['y'][ok]);M.append(a[ok])
    y=np.concatenate(Y);m=np.concatenate(M)
    return dict(**prf(y,(m>0).astype(int)),auroc=round(float(roc_auc_score(y,m)),3))
OUT={}
for name,Wn,epochs in (('C900',900,3),('C300',300,3),('C300s',300,1)):
    W,S=Wn,Wn//2
    Xtr,Ytr,_=windows(train); Xt=torch.tensor(Xtr,dtype=torch.long); Yt=torch.tensor(Ytr,dtype=torch.float32)
    # lr (B2) not used here; predict_record needs lr and kfeat on windows of this W
    lr=LogisticRegression(max_iter=500,C=1.0).fit(kfeat(Xtr),(Ytr.mean(1)>0.5).astype(int))
    OUT[name]=dict(W=Wn,epochs=epochs,n_windows=int(len(Xt)),seeds={})
    for seed in (0,1,2):
        torch.manual_seed(seed); np.random.seed(seed)
        net=Net(); opt=torch.optim.Adam(net.parameters(),1e-3); lossf=nn.BCEWithLogitsLoss(); steps=0
        for ep in range(epochs):
            perm=torch.randperm(len(Xt))
            for i in range(0,len(Xt),32):
                b=perm[i:i+32]; opt.zero_grad(); lossf(net(Xt[b]),Yt[b]).backward(); opt.step(); steps+=1
        net.eval()
        OUT[name]['seeds'][seed]=dict(steps=steps,cold=metr(test_cold),random=metr(test_rand)); print(name,seed,OUT[name]['seeds'][seed],flush=True)
        json.dump(OUT,open('results/ext_d6.json','w'),indent=1)
for name,v in OUT.items():
    for sp in ('cold','random'):
        ms=[v['seeds'][s][sp]['mcc'] for s in v['seeds']]; v[sp+'_mcc_mean']=round(float(np.mean(ms)),3); v[sp+'_mcc_sd']=round(float(np.std(ms)),3)
d=OUT['C300']['cold_mcc_mean']-OUT['C900']['cold_mcc_mean']; ds=OUT['C300s']['cold_mcc_mean']-OUT['C900']['cold_mcc_mean']
OUT['gate']=dict(delta_C300=round(d,3),delta_C300s=round(ds,3),met=bool(abs(d)>=0.05 and abs(ds)>=0.05 and d*ds>0),confounded_by_steps=bool(abs(d)>=0.05 and not(abs(ds)>=0.05 and d*ds>0)))
json.dump(OUT,open('results/ext_d6.json','w'),indent=1); print('done',OUT['gate'],flush=True)
