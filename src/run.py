import gzip,json,random,numpy as np,time,sys
from Bio import SeqIO
import torch,torch.nn as nn
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import f1_score,matthews_corrcoef,roc_auc_score
torch.set_num_threads(2); random.seed(0); np.random.seed(0); torch.manual_seed(0)
IDX={'A':0,'C':1,'G':2,'T':3}
recs=[]
try:
    for r in SeqIO.parse(gzip.open('/tmp/vir/part.gz','rt'),'genbank'):
        L=len(r.seq)
        if not (1000<=L<=20000): continue
        cds=[f for f in r.features if f.type=='CDS']
        if not cds: continue
        m=np.zeros(L,dtype=np.int8)
        for f in cds:
            for p in f.location.parts: m[int(p.start):int(p.end)]=1
        s=str(r.seq).upper(); x=np.array([IDX.get(c,-1) for c in s],dtype=np.int8)
        tax=r.annotations.get('taxonomy',[]); grp=tax[4] if len(tax)>4 else (tax[3] if len(tax)>3 else 'NA')
        ov=any(int(a.location.start)<int(b.location.end) and int(b.location.start)<int(a.location.end) for i,a in enumerate(cds) for b in cds[i+1:])
        recs.append(dict(id=r.id,x=x,y=m,grp=grp,L=L,ov=ov,gc=float(np.mean([c in 'GC' for c in s])),ncds=len(cds)))
except EOFError: pass
print('records',len(recs),flush=True)
R={'d1':dict(n_records=len(recs),n_groups=len(set(r['grp'] for r in recs)),mean_len=round(float(np.mean([r['L'] for r in recs]))),coding_fraction=round(float(np.mean(np.concatenate([r['y'] for r in recs]))),3),overlapping_fraction=round(float(np.mean([r['ov'] for r in recs])),3))}
grps=sorted(set(r['grp'] for r in recs)); rng=random.Random(0); rng.shuffle(grps); held=set(grps[:max(1,len(grps)//4)])
train_pool=[r for r in recs if r['grp'] not in held]; cold=[r for r in recs if r['grp'] in held]
rng.shuffle(train_pool); rng.shuffle(cold)
train=train_pool[:400]; test_cold=cold[:200]
rest=[r for r in train_pool[400:]]; test_rand=rest[:200]
R['d1'].update(n_train=len(train),n_test_cold=len(test_cold),n_test_random=len(test_rand),held_out_groups=len(held))
def orf_mask(x,Lmin):
    L=len(x); m=np.zeros(L,dtype=np.int8)
    for strand in (0,1):
        s=x if strand==0 else (3-x[::-1]); s=np.where(x[0:0].size or True,s,s)
        s=np.where(s<0,-1,s)
        for f in range(3):
            start=f; 
            codons=[(i,tuple(s[i:i+3])) for i in range(f,L-2,3)]
            cur=f
            for i,c in codons:
                if c in ((3,0,0),(3,2,0),(3,0,2)):  # TAA TAG TGA
                    if i-cur>=Lmin:
                        a,b=cur,i+3
                        if strand==1: a,b=L-b,L-a
                        m[a:b]=1
                    cur=i+3
            if L-cur>=Lmin and L-cur>0:
                a,b=cur,L
                if strand==1: a,b=L-b,L-a
                m[a:b]=1
    return m
def prf(y,p): return dict(f1=round(float(f1_score(y,p)),3),mcc=round(float(matthews_corrcoef(y,p)),3))
# baseline B1: pick L on train (subset for speed)
sub=train[:80]; best=None; R['d7']={}
for Lm in (150,300,450,600,900):
    yy=np.concatenate([r['y'] for r in sub]); pp=np.concatenate([orf_mask(r['x'],Lm) for r in sub]); f=f1_score(yy,pp); R['d7'][Lm]=round(float(f),3)
    if best is None or f>best[0]: best=(f,Lm)
Lbest=best[1]; R['b1_L_chosen_on_train']=Lbest; print('B1 L',Lbest,flush=True)
# windows
W,S=900,450
def windows(rs):
    X=[];Y=[];own=[]
    for k,r in enumerate(rs):
        L=len(r['x']); starts=list(range(0,max(1,L-W+1),S)); 
        if starts[-1]+W<L: starts.append(L-W)
        for a in starts:
            x=r['x'][a:a+W]; y=r['y'][a:a+W]
            if len(x)<W: continue
            X.append(np.where(x<0,0,x)); Y.append(y); own.append((k,a))
    return np.array(X),np.array(Y),own
Xtr,Ytr,_=windows(train); print('train windows',len(Xtr),flush=True)
# B2 kmer LR
def kfeat(X):
    K=np.zeros((len(X),64+1),dtype=np.float32); 
    c=X[:,:-2]*16+X[:,1:-1]*4+X[:,2:]
    for i in range(len(X)): K[i,:64]=np.bincount(c[i],minlength=64)/c.shape[1]; K[i,64]=np.isin(X[i],(1,2)).mean()
    return K
lr=LogisticRegression(max_iter=500,C=1.0).fit(kfeat(Xtr),(Ytr.mean(1)>0.5).astype(int))
class Net(nn.Module):
    def __init__(s):
        super().__init__(); s.e=nn.Embedding(4,16); s.c=nn.Conv1d(16,64,9,padding=4); s.g=nn.GRU(64,64,batch_first=True,bidirectional=True); s.o=nn.Linear(128,1)
    def forward(s,x):
        h=torch.relu(s.c(s.e(x).transpose(1,2))).transpose(1,2); h,_=s.g(h); return s.o(h).squeeze(-1)
net=Net(); opt=torch.optim.Adam(net.parameters(),1e-3); lossf=nn.BCEWithLogitsLoss()
Xt=torch.tensor(Xtr,dtype=torch.long); Yt=torch.tensor(Ytr,dtype=torch.float32); t0=time.time()
for ep in range(3):
    perm=torch.randperm(len(Xt)); tot=0
    for i in range(0,len(Xt),32):
        b=perm[i:i+32]; opt.zero_grad(); l=lossf(net(Xt[b]),Yt[b]); l.backward(); opt.step(); tot+=l.item()*len(b)
    print('epoch',ep,round(tot/len(Xt),4),round(time.time()-t0),'s',flush=True)
def predict_record(r):
    x=r['x']; L=len(x); X,_,own=windows([r]); 
    if len(X)==0: return None
    with torch.no_grad(): lo=net(torch.tensor(X,dtype=torch.long)).numpy()
    acc=np.zeros(L); cnt=np.zeros(L); 
    for (k,a),p in zip(own,lo): acc[a:a+W]+=p; cnt[a:a+W]+=1
    ok=cnt>0; acc[ok]/=cnt[ok]
    kf=lr.predict_proba(kfeat(X))[:,1]; b2=np.zeros(L); c2=np.zeros(L)
    for (k,a),p in zip(own,kf): b2[a:a+W]+=p; c2[a:a+W]+=1
    ok2=c2>0; b2[ok2]/=c2[ok2]; return acc,b2,ok
def evaluate(rs,name):
    Y=[];M=[];B=[];O=[];S=[];per=[]
    for r in rs:
        out=predict_record(r)
        if out is None: continue
        a,b,ok=out; o=orf_mask(r['x'],Lbest); y=r['y'][ok]
        Y.append(y);M.append(a[ok]);B.append(b[ok]);O.append(o[ok]); per.append(dict(grp=r['grp'],L=r['L'],gc=r['gc'],ov=r['ov'],f1_m=f1_score(y,(a[ok]>0).astype(int)) if y.sum()>0 else np.nan,f1_o=f1_score(y,o[ok]) if y.sum()>0 else np.nan,f1_b=f1_score(y,(b[ok]>0.5).astype(int)) if y.sum()>0 else np.nan, edge=None,
            err_edge=float(np.mean([min(np.abs(np.where(np.diff(r['y'][ok].astype(int))!=0)[0]-i).min() if (np.diff(r['y'][ok].astype(int))!=0).any() else 1e9,1e9)<=30 for i in np.where(((a[ok]>0).astype(int))!=y)[0][:300]])) if ((a[ok]>0).astype(int)!=y).any() else np.nan))
    y=np.concatenate(Y); m=np.concatenate(M); b=np.concatenate(B); o=np.concatenate(O)
    res=dict(n_records=len(per),n_bases=int(len(y)),coding_fraction=round(float(y.mean()),3),
      ORF=prf(y,o),KmerLR=dict(**prf(y,(b>0.5).astype(int)),auroc=round(float(roc_auc_score(y,b)),3)),Net=dict(**prf(y,(m>0).astype(int)),auroc=round(float(roc_auc_score(y,m)),3)))
    return res,per
for name,rs in (('d2_cold',test_cold),('d3_random',test_rand)):
    res,per=evaluate(rs,name); R[name]=res; R[name+'_per']=per if name=='d2_cold' else None; print(name,res,flush=True)
per=R['d2_cold_per']; import pandas as pd; P=pd.DataFrame(per); R['d2_cold_per']=None
best_base=max(R['d2_cold']['ORF']['f1'],R['d2_cold']['KmerLR']['f1']); R['d2_gate_net_beats_best_baseline_by_0.02']=bool(R['d2_cold']['Net']['f1']>=best_base+0.02)
R['d3_gate_gap_le_0.05']=bool(R['d3_random']['Net']['f1']-R['d2_cold']['Net']['f1']<=0.05)
g=P.groupby('grp').agg(n=('L','size'),f1_net=('f1_m','mean'),f1_orf=('f1_o','mean')).round(3); R['d4_per_group']=g.to_dict('index')
P['gc_bin']=pd.cut(P.gc,[0,0.4,0.5,0.6,1]); R['d5_gc']={str(k):dict(n=int(len(v)),f1_net=round(float(v.f1_m.mean()),3),f1_orf=round(float(v.f1_o.mean()),3)) for k,v in P.groupby('gc_bin',observed=True)}
R['d8_fraction_errors_within_30nt_of_cds_edge']=round(float(P.err_edge.mean()),3)
R['d9_overlap']={str(k):dict(n=int(len(v)),f1_net=round(float(v.f1_m.mean()),3),f1_orf=round(float(v.f1_o.mean()),3)) for k,v in P.groupby('ov')}
P['len_bin']=pd.cut(P.L,[0,3000,7000,12000,20001]); R['d10_length']={str(k):dict(n=int(len(v)),f1_net=round(float(v.f1_m.mean()),3),f1_orf=round(float(v.f1_o.mean()),3)) for k,v in P.groupby('len_bin',observed=True)}
json.dump(R,open('results/results.json','w'),indent=1,default=str); print(json.dumps(R,indent=1,default=str)[:5000])
