import json
S=open('src/run.py').read()
pre=S.split("net=Net(); opt")[0]
exec(pre)
mid=S.split("def predict_record")[1].split("for name,rs in (('d2_cold'")[0]
exec("def predict_record"+mid)
import pyrodigal
gf=pyrodigal.GeneFinder(meta=True)
def pyr_mask(r):
    L=len(r['x']); m=np.zeros(L,dtype=np.int8)
    seq=''.join('ACGT'[v] if v>=0 else 'N' for v in r['x'])
    for g in gf.find_genes(seq.encode()): m[g.begin-1:g.end]=1
    return m
def train(seed):
    global net
    torch.manual_seed(seed); np.random.seed(seed)
    net=Net(); opt=torch.optim.Adam(net.parameters(),1e-3); lossf=nn.BCEWithLogitsLoss()
    for ep in range(3):
        perm=torch.randperm(len(Xt))
        for i in range(0,len(Xt),32):
            b=perm[i:i+32]; opt.zero_grad(); lossf(net(Xt[b]),Yt[b]).backward(); opt.step()
    net.eval(); return net
Xt=torch.tensor(Xtr,dtype=torch.long); Yt=torch.tensor(Ytr,dtype=torch.float32)
def collect(rs):
    Y=[];M=[];O=[];Pm=[];G=[]
    for r in rs:
        out=predict_record(r)
        if out is None: continue
        a,b,ok=out; Y.append(r['y'][ok]);M.append(a[ok]);O.append(orf_mask(r['x'],Lbest)[ok]);Pm.append(pyr_mask(r)[ok]);G.append(np.full(ok.sum(),r['grp']))
    return [np.concatenate(z) for z in (Y,M,O,Pm,G)]
R={'Lbest':Lbest}
def m_(y,p): return dict(f1=float(f1_score(y,p)),mcc=float(matthews_corrcoef(y,p)))
R['A2']={}
for seed in (0,1,2):
    train(seed); res={}
    for nm,rs in (('cold',test_cold),('random',test_rand)):
        y,m,o,p,gg=collect(rs); res[nm]=dict(net=dict(**m_(y,(m>0).astype(int)),auroc=float(roc_auc_score(y,m))))
        if seed==0:
            res[nm]['orf']=m_(y,o); res[nm]['pyrodigal']=m_(y,p)
            if nm=='cold':
                k=gg!='Polyploviricotina'; res['cold_excl_Polyploviricotina']=dict(n_bases=int(k.sum()),net=m_(y[k],(m[k]>0).astype(int)),orf=m_(y[k],o[k]),pyrodigal=m_(y[k],p[k]),n_groups_in_cold=int(len(set(gg))),poly_fraction_of_bases=float(1-k.mean()))
    R['A2'][seed]=res; print('seed',seed,json.dumps(res),flush=True)
for nm in ('cold','random'):
    for k in ('f1','mcc','auroc'):
        v=[R['A2'][s][nm]['net'][k] for s in (0,1,2)]; R.setdefault('A2_summary',{}).setdefault(nm,{})[k]=dict(mean=float(np.mean(v)),sd=float(np.std(v)),values=v)
R['A1']={nm:dict(pyrodigal=R['A2'][0][nm]['pyrodigal'],orf=R['A2'][0][nm]['orf'],net_seed0=R['A2'][0][nm]['net']) for nm in ('cold','random')}
R['A1']['statement_net_loses_to_pyrodigal_cold']=bool(R['A1']['cold']['pyrodigal']['mcc']>R['A2_summary']['cold']['mcc']['mean']+0.05)
json.dump(R,open('results/audit.json','w'),indent=1); print('done',json.dumps(R['A1']))
