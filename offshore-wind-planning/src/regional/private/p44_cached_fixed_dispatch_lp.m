function [model,ix] = p44_cached_fixed_dispatch_lp(T,sp)
% Original B0 sparse LP; callers bind Emax to the installed capacity.
persistent cachedSignature cachedModel cachedIndex
signature=[T,sp.eta_ch,sp.eta_dis,sp.self_discharge_per_h, ...
    sp.soc_min,sp.soc_max,sp.duration_h];
if isequal(signature,cachedSignature)
    model=cachedModel; ix=cachedIndex; return;
end
assert(sp.duration_h>0&&sp.eta_ch>0&&sp.eta_ch<=1&&sp.eta_dis>0&&sp.eta_dis<=1);
assert(sp.soc_min>=0&&sp.soc_min<sp.soc_max&&sp.soc_max<=1);
assert(sp.self_discharge_per_h>=0&&sp.self_discharge_per_h<1);
ix.Eini=1; ix.Emax=2;
ix.charge=(3:T+2).'; ix.discharge=(T+3:2*T+2).';
ix.energy=(2*T+3:3*T+3).'; ix.delivery=(3*T+4:4*T+3).';
ix.curtail=(4*T+4:5*T+3).'; nv=5*T+3; nr=6*T+6;
rr=[]; cc=[]; vv=[];
    function add(ar,ac,av)
        rr=[rr;ar(:)]; cc=[cc;ac(:)]; vv=[vv;av(:)];
    end
q=(1:T).'; u=ones(T,1);
add(q,ix.charge,u); add(q,repmat(ix.Emax,T,1),-u/sp.duration_h);
add(T+q,ix.discharge,u); add(T+q,repmat(ix.Emax,T,1),-u/sp.duration_h);
add(2*T+1,ix.energy(1),1); add(2*T+1,ix.Eini,-1);
add(2*T+2,ix.Emax,sp.soc_min); add(2*T+2,ix.Eini,-1);
add(2*T+3,ix.Eini,1); add(2*T+3,ix.Emax,-sp.soc_max);
ix.balanceRows=2*T+3+q;
add(ix.balanceRows,ix.discharge,u); add(ix.balanceRows,ix.charge,-u);
add(ix.balanceRows,ix.delivery,-u); add(ix.balanceRows,ix.curtail,-u);
recur=3*T+3+q;
add(recur,ix.energy(2:end),u); add(recur,ix.energy(1:end-1),-(1-sp.self_discharge_per_h)*u);
add(recur,ix.charge,-sp.eta_ch*u); add(recur,ix.discharge,u/sp.eta_dis);
qs=(1:T+1).'; us=ones(T+1,1); lower=4*T+3+qs; upper=5*T+4+qs;
add(lower,repmat(ix.Emax,T+1,1),sp.soc_min*us); add(lower,ix.energy,-us);
add(upper,ix.energy,us); add(upper,repmat(ix.Emax,T+1,1),-sp.soc_max*us);
add(nr,ix.Eini,1); add(nr,ix.energy(end),-1);
model=struct('A',sparse(rr,cc,vv,nr,nv),'rhs',zeros(nr,1), ...
    'sense',repmat('<',nr,1),'obj',zeros(nv,1),'lb',zeros(nv,1), ...
    'ub',inf(nv,1),'vtype',repmat('C',nv,1),'modelsense','min', ...
    'modelname','P44_fixed_original_B0_dispatch');
model.sense(2*T+1)='='; model.sense(ix.balanceRows)='='; model.sense(recur)='=';
cachedSignature=signature; cachedModel=model; cachedIndex=ix;
end
