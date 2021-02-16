% program for analysis of Tasmanian devil population projection

% ----------------------------------------
% model with Guiler (1978) survival rates
% ----------------------------------------

format long g

% define max age
age_max=6;

% define survival
s0g=0.3983; % Guiler (1978)
s1g=0.6278; % Guiler (1978)
s2g=0.6524; % Guiler (1978)
s3g=0.6220; % Guiler (1978)
s4g=0.5285; % Guiler (1978)
s5g=0.2711; % Guiler (1978)

% define fecundity parameters
m_imm=0.061; % 1 female seen with pouch young (Pemberton 1990)
m_primi=2.0; % average pouch young per primiparous female (Guiler 1978)
m_prime=3.6; % average pouch young per 'prime' (3-5) female (Guiler 1978)
m_old=2.0; % average pouch young per old (6) female (Guiler 1978)

rs_lo=0.80; % proportion of females carrying young (Guiler 1978)
%rs_hi=0.81; % proportion of females carrying young (Pemberton 1990)
rs=rs_lo;

primi=2; % age at first reproduction

% initial population sizes
Nmin=130000; % N. Mooney, pers. comm (2003)
Nmax=150000; % N. Mooney, pers. comm (2003) 
Navgmax = mean([Nmin Nmax]);
Navg50 = Navgmax / 2; % start out 1/2 population size

% sex ratios
sr_tot1 = 0.5 % proportion female - Pemberton (1990)
sr_tot2 = 0.5362 % Guiler (1970a)
sr_avg = mean([sr_tot1,sr_tot2]);

x1=0.48; % pouch young proportion female (Guiler 1978)
x2=0.58; % pouch young proportion female (Guiler 1970a)
x3=0.57; % pouch young proportion female (Guiler 1970a)
x4=0.56; % pouch young proportion female (Hughes 1982)
x5=0.50; % pouch young proportion female (Pemberton 1990)
x_avg = mean([x1,x2,x3,x4,x5]);

sr = sr_tot2
x=x_avg

%total females in population
f=Navg50/2;

%Initial population size vector
N=f;

% the normal matrix
ag=[0 s0g*m_primi*x*rs s0g*m_prime*x*rs s0g*m_prime*x*rs s0g*m_prime*x*rs s0g*m_old*x*rs
    s1g 0 0 0 0 0
    0 s2g 0 0 0 0
    0 0 s3g 0 0 0
    0 0 0 s4g 0 0
    0 0 0 0 s5g 0];

% eigenvalues & eigenvectors
[w,d] = eig(ag);
v=conj(inv(w));

%lambda
lambda=diag(d);

% max lambdas
[maxlambda,imax]=max(diag(d));

%logmaxlambda=exponential rate of increase (r)
r=log(maxlambda);

% stable age distribution
w=w(:,imax);

% reproductive values
v=real(v(imax,:))';

% sensitivities and elasticities
senmat=(v*w')/(v'*w);
emat=(ag.*senmat)/maxlambda;

%damping ratios
rho=lambda(1)./abs(lambda);

% periods of oscillation
period=2*pi./angle(lambda);

% stable size distribution (scaled to sum to 1)
ssd_g=w(:,1)/sum(w(:,1));

% ssd classes
ssd_juv = sum(ssd_g(1:1,1));
ssd_ad = sum(ssd_g(2:6,1));

% reproductive value (scaled so stage 1 = 1)
reprovalue=real(v(1,:)')/real(v(1,1));

% mean generation time (avg age breeding females)
%gen=log(0.9846)/(log(maxlambda));

% size of matrix
k=size(ag,1);

% pick initial vectors
ng=ones(k,1);
%n(1:2,1)=0.2382*N/2; % Guiler 1978
%n(3:6,1)=0.7618*N/4; % Guiler 1978
ng=ssd_g*N;

% age vector
age_vec=zeros(k,1);
for a=1:k-1
    age_vec(a+1,1)=a;
end % a loop
age_vec=age_vec+1;

%Calculate Quasi-extinction times
Qg=(log(50/sum(ng)))/log(maxlambda);
if Qg < 0;
    Qg='infinity';
else Qg=Qg;
end

% do a simulation
% first specify the initial condition and length of simulation
tlimit=25;

%set population size year step vector
pop_vecg=ones(tlimit+1,1);
pop_vecg(1,1)=sum(ng);

%set year step vector
yr_vec=ones(tlimit+1,1);
for c=1:tlimit
    yr_vec(c+1,1)=c;
    yr_vec(1,1)=0;
end % c loop

% survival variance vectors
s1_vec = ones(tlimit+1,1);
s1_vec(1,1) = (s1g*(1-s1g))/ng(1,1);

s2_vec = ones(tlimit+1,1);
s2_vec(1,1) = (s2g*(1-s2g))/ng(1,1);

% assume stable-age distribution
N_stable = ssd_g*N;

%then iterate
for i=1:tlimit;
    ng=ag*ng;
    s1_vec(i+1,1)=(s1g*(1-s1g))/ng(1,1);
    s2_vec(i+1,1)=(s2g*(1-s2g))/ng(2,1);
    pop_vecg(i+1,1)=(sum(ng));    
end

log_pop_vecg=log10(pop_vecg);

%total population size after 'tlimit' years
pop_st=N
pop_end=sum(ng)

%total population size after 'tlimit' years
N
tlimit
pop_end=sum(ng)
maxlambda
r
Qg
format short

% construct survival vector for display
survg=ones(k,1);
survg(1,1)=s0g;
survg(2,1)=s1g;
survg(3,1)=s2g;
survg(4,1)=s3g;
survg(5,1)=s4g;
survg(6,1)=s5g;

% construct fertility vector for display
fertg=ones(k,1);
fertg(1,1)=m_imm*x;
fertg(2,1)=m_primi*x*rs;
fertg(3,1)=m_prime*x*rs;
fertg(4,1)=m_prime*x*rs;
fertg(5,1)=m_prime*x*rs;
fertg(6,1)=m_old&x*rs;

% continue displays
ag
ssd_juv
ssd_ad
emat

%Make density independent plots
subplot(2,2,1), plot(yr_vec,pop_vecg);
axis square;
axis([-0.5 tlimit+1 0 (max(pop_vecg)+(0.25*(max(pop_vecg))))]);
xlabel('year');
ylabel('N females');
subplot(2,2,2), plot(yr_vec,log_pop_vecg);
axis square;
axis([-0.5 tlimit+1 0 (max(log_pop_vecg)+(0.25*(max(log_pop_vecg))))]);
xlabel('year');
ylabel('log N females');

%Make survival and fecundity plots
subplot(2,2,3), plot(age_vec,survg);
axis square;
axis([0.5 k+0.5 0 1]);
xlabel('age in years');
ylabel('survival probability');
subplot(2,2,4), plot(age_vec,fertg);
axis square;
axis([0.5 k+0.5 0 2]);
xlabel('age in years');
ylabel('m (fertility)');

pause

% -----------------------------------------------------------
% model with U-shaped mortality to produce 'max' growth rate
% -----------------------------------------------------------

format long g

% define max age
age_max=6;

% re-define survival

s0_hi=0.500; % Pemberton (1990) - 1983/84 (weanlings)
s0_lo=0.360; % Pemberton (1990) - 1984/85 (weanlings)
s0_avg = mean([s0_hi,s0_lo])
s0 = s0_avg;

s1_lo=0.170; % Pemberton (1990) - mean annual rates - lower
s1_hi=0.300;

s2_lo=0.820; % Pemberton (1990) - mean annual rates - lower
s2_hi=0.90;

%s1=mean([s1_lo,s1_hi]);
%s2=mean([s2_lo,s2_hi]);

s2=0.82 % Pemberton (1990), Method 3

s1=s1g; % using Guiler's data for U-shaped curve
s3=s2; % assuming no change in survival until last 1 years
s4=s2; % assuming no change in survival until last 1 years
s5=s5g; % using Guiler's data for U-shaped curve

% define fecundity parameters
m_imm=0.061; % 1 female seen with pouch young (Pemberton 1990)
m_primi=2.0; % average pouch young per primiparous female (Guiler 1978)
m_prime=3.6; % average pouch young per 'prime' (3-5) female (Guiler 1978)
m_old=2.0; % average pouch young per old (6) female (Guiler 1978)

%rs_lo=0.80; % proportion of females carrying young (Guiler 1978)
rs_hi=0.81; % proportion of females carrying young (Pemberton 1990)
rs=rs_hi;

primi=2; % age at first reproduction

% initial population sizes
Nmin=130000; % N. Mooney, pers. comm (2003)
Nmax=150000; % N. Mooney, pers. comm (2003) 
Navgmax = mean([Nmin Nmax]);
Navg50 = Navgmax / 2; % start out 1/2 population size

% sex ratios
sr_tot1 = 0.5 % proportion female - Pemberton (1990)
sr_tot2 = 0.5362 % Guiler (1970a)
sr_avg = mean([sr_tot1,sr_tot2]);

x1=0.48; % pouch young proportion female (Guiler 1978)
x2=0.58; % pouch young proportion female (Guiler 1970a)
x3=0.57; % pouch young proportion female (Guiler 1970a)
x4=0.56; % pouch young proportion female (Hughes 1982)
x5=0.50; % pouch young proportion female (Pemberton 1990)
x_avg = mean([x1,x2,x3,x4,x5]);

sr = sr_avg
x=x_avg

%total females in population
f=Navg50/2;

%Initial population size vector
N=f;

% the normal matrix
amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
       s1 0 0 0 0 0
       0 s2 0 0 0 0
       0 0 s3 0 0 0
       0 0 0 s4 0 0
       0 0 0 0 s5 0];
   
% eigenvalues & eigenvectors
[w,d] = eig(amax);
v=conj(inv(w));

%lambda
lambda=diag(d);

% max lambdas
[maxlambda,imax]=max(diag(d));

%logmaxlambda=exponential rate of increase (r)
r=log(maxlambda);

% stable age distribution
w=w(:,imax);

% reproductive values
v=real(v(imax,:))';

% sensitivities and elasticities
senmat=(v*w')/(v'*w);
emat=(amax.*senmat)/maxlambda;

%damping ratios
rho=lambda(1)./abs(lambda);

% periods of oscillation
period=2*pi./angle(lambda);

% stable size distribution (scaled to sum to 1)
ssd_inc=w(:,1)/sum(w(:,1));

% ssd classes
ssd_juv = sum(ssd_inc(1:1,1));
ssd_ad = sum(ssd_inc(2:6,1));

% reproductive value (scaled so stage 1 = 1)
reprovalue=real(v(1,:)')/real(v(1,1));

% mean generation time (avg age breeding females)
%gen=log(0.9846)/(log(maxlambda));

% size of matrix
k=size(amax,1);

% pick initial vectors
n=ones(k,1);
%n(1:2,1)=0.2382*N/2; % Guiler 1978
%n(3:6,1)=0.7618*N/4; % Guiler 1978
n=ssd_inc*N;

% age vector
age_vec=zeros(k,1);
for a=1:k-1
    age_vec(a+1,1)=a;
end % a loop
age_vec=age_vec+1;

%Calculate Quasi-extinction times
Q=(log(50/sum(n)))/log(maxlambda);
if Q < 0;
    Q='infinity';
else Q=Q;
end

% do a simulation
% first specify the initial condition and length of simulation
tlimit=25;

%set population size year step vector
pop_vec=ones(tlimit+1,1);
pop_vec(1,1)=sum(n);

%set year step vector
yr_vec=ones(tlimit+1,1);
for c=1:tlimit
    yr_vec(c+1,1)=c;
    yr_vec(1,1)=0;
end % c loop

% survival variance vectors
s1_vec = ones(tlimit+1,1);
s1_vec(1,1) = (s1*(1-s1))/n(1,1);

s2_vec = ones(tlimit+1,1);
s2_vec(1,1) = (s2*(1-s2))/n(1,1);

% assume stable-age distribution
N_stable = ssd_inc*N;

%then iterate
for i=1:tlimit;
    n=amax*n;
    s1_vec(i+1,1)=(s1*(1-s1))/n(1,1);
    s2_vec(i+1,1)=(s2*(1-s2))/n(2,1);
    pop_vec(i+1,1)=(sum(n));    
end

log_pop_vec=log10(pop_vec);

%total population size after 'tlimit' years
pop_st=N
pop_end=sum(n)

%total population size after 'tlimit' years
N
tlimit
pop_end=sum(n)
maxlambda
r
Q
format short

% construct survival vector for display
surv=ones(k,1);
surv(1,1)=s0;
surv(2,1)=s1;
surv(3,1)=s2;
surv(4,1)=s3;
surv(5,1)=s4;
surv(6,1)=s5;

% construct fertility vector for display
fert=ones(k,1);
fert(1,1)=m_imm*x;
fert(2,1)=m_primi*x*rs;
fert(3,1)=m_prime*x*rs;
fert(4,1)=m_prime*x*rs;
fert(5,1)=m_prime*x*rs;
fert(6,1)=m_old&x*rs;

% continue displays
amax
ssd_juv
ssd_ad
emat

%Make density independent plots
subplot(2,2,1), plot(yr_vec,pop_vec);
axis square;
axis([-0.5 tlimit+1 0 (max(pop_vec)+(0.25*(max(pop_vec))))]);
xlabel('year');
ylabel('N females');
subplot(2,2,2), plot(yr_vec,log_pop_vec);
axis square;
axis([-0.5 tlimit+1 0 (max(log_pop_vec)+(0.25*(max(log_pop_vec))))]);
xlabel('year');
ylabel('log N females');

%Make survival and fecundity plots
subplot(2,2,3), plot(age_vec,surv);
axis square;
axis([0.5 k+0.5 0 1]);
xlabel('age in years');
ylabel('survival probability');
subplot(2,2,4), plot(age_vec,fert);
axis square;
axis([0.5 k+0.5 0 2]);
xlabel('age (years)');
ylabel('m (fertility)');


% ----------------------------------------------------------------------------
% deterministic model with negative density feedback in survival (no disease)
% ----------------------------------------------------------------------------

% Re-define new N females
Nnew = 35000;

% set new time limit
tlimit = 200;

% reset amax and adisease to original values
amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
       s1 0 0 0 0 0
       0 s2 0 0 0 0
       0 0 s3 0 0 0
       0 0 0 s4 0 0
       0 0 0 0 s5 0];

addd=amax;

    % pick initial vectors
    nddd=ones(k,1);
    nddd=ssd_inc*Nnew;
    
    %set year step vector
    yr_vecddd=ones(tlimit+1,1);
    for c=1:tlimit
        yr_vecddd(c+1,1)=c;
        yr_vecddd(1,1)=0;
    end % c loop

    %set population size year step vector
    pop_vecddd=ones(tlimit+1,1);
    pop_vecddd(1,1)=sum(nddd);

    % set dd survival vector
    sddd_vec = zeros(tlimit+1,1);
    sddd_vec(1,1)=s2;
    
    %then iterate
    for iy=1:tlimit;

        nddd=addd*nddd;
    
        % add to population vector
        pop_vecddd(iy+1,1)=(sum(nddd));
    
        % Set negative density feedback function (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))

        % parameters with 35000 mid-point
        acoeffd=0.1968;
        bcoeffd=2.9838;
        x0d=70000;
        y0d=0.6234;

        % Redefine survival probabilities
        s2_ddd = y0d+(acoeffd/(1+(sum(nddd)/x0d)^bcoeffd));
        s3_ddd = s2_ddd;
        s4_ddd = s2_ddd;
        s5_ddd = s2_ddd;
    
        % the new matrix
        addd=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
              s1 0 0 0 0 0
              0 s2_ddd 0 0 0 0
              0 0 s3_ddd 0 0 0
              0 0 0 s4_ddd 0 0
              0 0 0 0 s5_ddd 0];

        % add updated survival probability to survival vector
        sddd_vec(iy+1,1)=s2_ddd;
          
    end % for iy

    % calculate instantaneous growth rates
    lpvddd = length(pop_vecddd);
    rddd = pop_vecddd(2:lpvddd) ./ pop_vecddd(1:lpvddd-1);
    lrddd = log(rddd);
    
    % log-transform population vector
    log_pop_vecddd=log10(pop_vecddd);

    % number of years to stability
    diff_thresh = 10; % define 'stability'
    pop_diff = diff(pop_vecddd);
    [Iddd,Jddd] = find(abs(pop_diff)<diff_thresh);
    clear Jddd;
    stab_yrs = Iddd(1) % years to 'stability'

    pop_vecddd(tlimit+1) % final population size
    
    sddd_vec(tlimit+1) % final adult survival at stability
    
    % Plot trajectory
    subplot(1,3,1), plot(yr_vecddd,pop_vecddd);
    axis square;
    axis([-0.5 tlimit+1 0 (max(pop_vecddd)+(0.25*(max(pop_vecddd))))]);
    xlabel('year');
    ylabel('N females');
    hold
    
    % create stability line
    popmax = (max(pop_vecddd)+(0.25*(max(pop_vecddd))));
    plot(stab_yrs,0:popmax,'r-');    
    hold
    
    subplot(1,3,2), plot(yr_vecddd,sddd_vec);
    axis square;
    axis([-0.5 tlimit+1 0.5 1]);
    xlabel('year');
    ylabel('p adult survival');
    hold
    plot(stab_yrs,0:1,'r:');    
    hold
    
    subplot(1,3,3), plot(yr_vecddd(1:tlimit),lrddd);
    axis square;
    axis([-0.5 tlimit+1 min(lrddd)-(0.1*(max(lrddd))) max(lrddd)+0.1*(abs(min(lrddd))) ]);
    xlabel('year');
    ylabel('r');


% ----------------------------------------
% model with disease-related mortality
% ----------------------------------------

% lower survival with disease
s_ad_disease_lo = 0.10;
s_ad_disease_hi = 0.60;

s_ad_disease_avg = mean([s_ad_disease_lo,s_ad_disease_hi]);
%s_ad_disease_avg = 0.7

sc2=s_ad_disease_avg;
sc3=s_ad_disease_avg;
sc4=s_ad_disease_avg;
sc5=s_ad_disease_avg;

% the disease matrix
adisease=[0 s0*m_imm*x s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
         s1 0 0 0 0 0
         0 sc2 0 0 0 0
         0 0 sc3 0 0 0
         0 0 0 sc4 0 0
         0 0 0 0 sc5 0];

% eigenvalues & eigenvectors
[wc,dc] = eig(adisease);
vc=conj(inv(wc));

%lambda
lambdac=diag(dc);

% max lambdas
[maxlambdac,imax]=max(diag(dc));

%logmaxlambda=exponential rate of increase (r)
rc=log(maxlambdac);

% stable age distribution
wc=wc(:,imax);

% reproductive values
vc=real(vc(imax,:))';

% sensitivities and elasticities
senmatc=(vc*wc')/(vc'*wc);
ematc=(adisease.*senmatc)/maxlambdac;

% stable age distribution
w=w(:,imax);

% reproductive values
v=real(v(imax,:))';

%damping ratios
rho_c=lambda(1)./abs(lambdac);

% periods of oscillation
period_c=2*pi./angle(lambdac);

% stable size distribution (scaled to sum to 1)
ssd_c=wc(:,1)/sum(wc(:,1));

% ssd classes
ssdc_juv = sum(ssd_c(1:1,1));
ssdc_ad = sum(ssd_c(2:6,1));

% reproductive value (scaled so stage 1 = 1)
reprovalue=real(v(1,:)')/real(v(1,1));


% pick initial vectors
nc=ones(k,1);
%nc(1:2,1)=0.2382*N/2; % Guiler 1978
%nc(3:6,1)=0.7618*N/4; % Guiler 1978
nc=ssd_inc*N;

% half of population
N_half = N/2;

%Calculate Quasi-extinction times
Qc_half=(log(N_half/sum(nc)))/log(maxlambdac);
if Qc_half < 0;
    Qc_half='infinity';
else Qc_half=Qc_half;
end

% set disease tlimit
tlimitc = 6;

%set year step vector
yr_vecc=ones(tlimitc+1,1);
for c=1:tlimitc
    yr_vecc(c+1,1)=c;
    yr_vecc(1,1)=0;
end % c loop

%set population size year step vector
pop_vecc=ones(tlimitc+1,1);
pop_vecc(1,1)=sum(nc);

%then iterate
for i=1:tlimitc;
    nc=adisease*nc;
    pop_vecc(i+1,1)=(sum(nc));    
end

log_pop_vecc=log10(pop_vecc);

%total population size after 'tlimit' years
pop_stc=N
pop_endc=sum(nc)

%total population size after 'tlimit' years
N
tlimitc
pop_endc=sum(nc)
maxlambdac
rc
Qc_half
format short
adisease
ematc
nc
ssdc_juv
ssdc_ad

%Make density independent plots
subplot(1,2,1), plot(yr_vecc,pop_vecc);
axis square;
axis([-0.5 tlimitc+1 0 (max(pop_vecc)+(0.25*(max(pop_vecc))))]);
xlabel('year');
ylabel('N females (disease)');
hold
N_half_vec=ones(tlimitc+1,1)*N_half;
subplot(1,2,1), plot(yr_vecc,N_half_vec,'r:');
hold
subplot(1,2,2), plot(yr_vecc,log_pop_vecc);
axis square;
axis([-0.5 tlimitc+1 0 (max(log_pop_vecc)+(0.25*(max(log_pop_vecc))))]);
xlabel('year');
ylabel('log N females (disease)');

% ----------------------------------------------------------------------------------------------------------------------------------
% stochastic model to build disease bouts (Markhov Chain dependency) with neg dens feedback in survival - estimate p_dep & p_disease
% ----------------------------------------------------------------------------------------------------------------------------------

% reset amax and adisease to original values
amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
       s1 0 0 0 0 0
       0 s2 0 0 0 0
       0 0 s3 0 0 0
       0 0 0 s4 0 0
       0 0 0 0 s5 0];

adisease=[0 s0*m_imm*x s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
          s1 0 0 0 0 0
          0 sc2 0 0 0 0
          0 0 sc3 0 0 0
          0 0 0 sc4 0 0
          0 0 0 0 sc5 0];

adisdd = adisease;
      
% specify a time limit
tlimit=200;

% specify range of p_disease values
parraypdd = linspace(-1,1,40);

% specify range of p_dep values
parraydd=linspace(0,1,20);

% set storage matrices
all_loglambdasdd = zeros(40,20);
all_extdd = zeros(40,20);
all_minpcdd = zeros(40,20);

% set irpdd loop
for irpdd = 1:length(parraypdd);
    
    % set probability of p_disease
    p_disease = parraypdd(irpdd);


% set irdd loop
for irdd = 1:length(parraydd);

    % set probability of inter-year disease dependency
    p_depdd = parraydd(irdd);

    % Number of iterations to do
    iter = 1000;

% define matrix to hold population sizes
pop_vecdd_ci = ones(tlimit+1,iter);
log_vecdd_ci = ones(tlimit+1,iter);

% set extinction sum vector
  extinct=zeros(iter,1);

% do loop to estimate population size confidence intervals
for idd=1:iter;

    % generate random number tlimit vector
    random=rand(1,tlimit);

    % initial stochastic determination (1=disease; 2=no disease);
    ydd=(p_disease<=random);
    disdd=ydd+1;

    % start Markhov Chain
    disdd_new=disdd;

    for ip=1:tlimit-1;
    
        if (disdd_new(ip)==1) & (rand(1,1) <= p_depdd)
             disdd_new(ip+1) = 1;
        else
             disdd_new(ip+1) = disdd_new(ip+1);
        end

    end % for ip

    % pick initial vectors
    ndd=ones(k,1);
    %ndd(1:2,1)=0.2382*N/2; % Guiler 1978
    %ndd(3:6,1)=0.7618*N/4; % Guiler 1978
    ndd=ssd_inc*N;
    
    %set year step vector
    yr_vecdd=ones(tlimit+1,1);
    for c=1:tlimit
        yr_vecdd(c+1,1)=c;
        yr_vecdd(1,1)=0;
    end % c loop

    %set population size year step vector
    pop_vecdd=ones(tlimit+1,1);
    pop_vecdd(1,1)=sum(ndd);

    % set population extinction indicator vector
    ext_vecdd=zeros(tlimit+1,1);
    ext_vecdd(1,1)=0;

    %then iterate
    for iy=1:tlimit;

        if (disdd_new(iy)==1)
            ndd=adisdd*ndd;
        else
            ndd=amax*ndd;
        end

        % add to population vector
        pop_vecdd(iy+1,1)=(sum(ndd));

        % has the population gone extinct (i.e., n vector sums to < 1)?
        ext_thresh = 50; % set extinction threshold

        if sum(ndd) < ext_thresh;
            ext_vecdd(iy+1,1)=1;
        else
            ext_vecdd(iy+1,1)=0;
        end

        % Set negative density feedback function for the 'normal' matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))

        % parameters with 70000 mid-point
        acoeffn=0.1968;
        bcoeffn=2.9838;
        x0n=70000;
        y0n=0.6234;

        % Redefine survival probabilities
        s2_n = y0n+(acoeffn/(1+(sum(ndd)/x0n)^bcoeffn));
        s3_n = s2_n;
        s4_n = s2_n;
        s5_n = s2_n;
    
        % the new matrix
        amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
              s1 0 0 0 0 0
              0 s2_n 0 0 0 0
              0 0 s3_n 0 0 0
              0 0 0 s4_n 0 0
              0 0 0 0 s5_n 0];

        % Set negative density feedback function for the disease matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))
        acoeffdd=0.5629868;
        bcoeffdd=3.00328;
        x0dd=35000;
        y0dd=0.037572;

        % Redefine survival probabilities
        sc1_dd = s1;
        sc2_dd = y0dd+(acoeffdd/(1+(sum(ndd)/x0dd)^bcoeffdd));
        sc3_dd = sc2_dd;
        sc4_dd = sc2_dd;
        sc5_dd = sc2_dd;
    
        % the new disease matrix
        adisdd=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
        sc1_dd 0 0 0 0 0
        0 sc2_dd 0 0 0 0
        0 0 sc3_dd 0 0 0
        0 0 0 sc4_dd 0 0
        0 0 0 0 sc5_dd 0];

    end % for iy

    % log-transform population vector
    log_pop_vecdd=log10(pop_vecdd);

    % place population size vectors into storage matrix
    pop_vecdd_ci(:,idd) = pop_vecdd;
    log_vecdd_ci(:,idd) = log_pop_vecdd;

   % did the population go extinct during this run?
    extinctdd(idd,1) = sum(ext_vecdd);

end % idd loop

% Calculate mean population sizes
mean_popdd = mean(pop_vecdd_ci,2);
mean_log_popdd = mean(log_vecdd_ci,2);

% calculate stochastic r
mean_popdd1=zeros(length(mean_popdd)+1,1);
mean_popdd1(2:length(mean_popdd)+1,1)=mean_popdd;
rdd = mean_popdd(2:length(mean_popdd),1) ./ mean_popdd1(2:length(mean_popdd),1);

% calculate lowest percentage of original population size
pcpop_vecdd_ci = (pop_vecdd_ci/N);
minpcdd = min(pcpop_vecdd_ci,[],1);
mean_minpcdd = (mean(minpcdd))*100;

% Calculate probability of extinction
extdd = zeros(iter,1);     
for ie = 1:iter;
        if extinctdd(ie,1) > 0;
            extdd(ie,1) = 1;
        else
            extdd(ie,1) = 0;
        end
    end

pr_extdd = (sum(extdd)) / iter;

% stochastic growth rate
    loglambdasdd(irdd) = mean(log(rdd));

% create probability of extinction vector
    pr_ext_vecdd(irdd) = pr_extdd;

% create mean minimum per cent of start N vector
    mean_minpc_vecdd(irdd) = mean_minpcdd;

irdd

end %irdd loop

% store values in respective rows in storage matrices
all_loglambdasdd(irpdd,:) = loglambdasdd;
all_extdd(irpdd,:) = pr_ext_vecdd;
all_minpcdd(irpdd,:) = mean_minpc_vecdd;

irpdd

end %irpdd loop

% Plot loglambdas as a function of p_dep & p_disease
mesh(all_loglambdasdd);
view(-123,15)
shading interp;
colormap([0 0 0]);
axis tight
set(gca,'XTickLabels',[0.25 0.50 0.75 1.00]);
set(gca,'YTickLabels',[-0.75 -0.50 -0.25 0.00 0.25 0.50 0.75 1.00]);
zlabel('log lambda');
xlabel('p disease occurrence');
ylabel('temporal autocorrelation');
hold
lab = zeros(40,20);
surf(lab);
shading interp;

% Find where loglambdas are ~ 0
nearz_low = -0.0005;
nearz_hi = 0.0005;
[I,J] = find(all_loglambdasdd > nearz_low & all_loglambdasdd < nearz_hi); % finds coordinates
p_disease_rangedd = parraypdd(I)
p_dep_rangedd = parraydd(J)
min_p_disease = min(p_disease_rangedd)
max_p_disease = max(p_disease_rangedd)
min_p_dep = min(p_dep_rangedd)
max_p_dep = max(p_dep_rangedd)

% Plot pr_ext as a function of p_dep & p_disease
mesh(all_extdd);
view(56,21)
shading flat;
colormap([0 0 0]);
axis tight
set(gca,'XTickLabels',[0.25 0.50 0.75 1.00]);
set(gca,'YTickLabels',[-0.75 -0.50 -0.25 0.00 0.25 0.50 0.75 1.00]);
zlabel('p quasi-extinction');
xlabel('p disease occurrence');
ylabel('p temporal autocorrelation');

% Plot min_pc as a function of p_dep & p_disease
mesh(all_minpcdd);
view(-107,12)
shading flat;
colormap([0 0 0]);
axis tight
set(gca,'XTickLabels',[0.25 0.50 0.75 1.00]);
set(gca,'YTickLabels',[-0.75 -0.50 -0.25 0.00 0.25 0.50 0.75 1.00]);
zlabel('mean minimum %');
xlabel('p disease occurrence');
ylabel('p temporal autocorrelation');

% -------------------------------------------------------------------------------------------------------------------------------------------
% stochastic model to build disease bouts (Markhov Chain dependency) with neg dens feedback - use p_dep & p_disease range from previous model
% -------------------------------------------------------------------------------------------------------------------------------------------

% reset amax and adisease to original values
amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
       s1 0 0 0 0 0
       0 s2 0 0 0 0
       0 0 s3 0 0 0
       0 0 0 s4 0 0
       0 0 0 0 s5 0];

adisease=[0 s0*m_imm*x s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
          s1 0 0 0 0 0
          0 sc2 0 0 0 0
          0 0 sc3 0 0 0
          0 0 0 sc4 0 0
          0 0 0 0 sc5 0];
      
adisdd=adisease;

% define range of probability of inter-year disease dependency
p_depdd_dif = max_p_dep - min_p_dep;

% define range of probability of p_disease
p_diseasedd_dif = max_p_disease - min_p_disease;

% define matrix to hold population sizes
pop_vecdd_ci = ones(tlimit+1,iter);
log_vecdd_ci = ones(tlimit+1,iter);

% set extinction sum vector
  extinct=zeros(iter,1);

% do loop to estimate population size confidence intervals
for idd=1:iter;

    % generate random number tlimit vector
    random=rand(1,tlimit);

    % sample p_disease from possible range (defined under previous model)
    p_diseasedd = (rand * p_diseasedd_dif) + min_p_disease;

    % initial stochastic determination (1=disease; 2=no disease);
    ydd=(p_diseasedd<=random);
    disdd=ydd+1;

   % sample p_dep from possible range (defined under previous model)
    p_depdd = (rand * p_depdd_dif) + min_p_dep;

    % start Markhov Chain
    disdd_new=disdd;

    for ip=1:tlimit-1;
    
        if (disdd_new(ip)==1) & (rand(1,1) <= p_depdd)
             disdd_new(ip+1) = 1;
        else
             disdd_new(ip+1) = disdd_new(ip+1);
        end

    end % for ip

    % pick initial vectors
    ndd=ones(k,1);
    %ndd(1:2,1)=0.2382*N/2; % Guiler 1978
    %ndd(3:6,1)=0.7618*N/4; % Guiler 1978
    ndd=ssd_inc*N;
    
    %set year step vector
    yr_vecdd=ones(tlimit+1,1);
    for c=1:tlimit
        yr_vecdd(c+1,1)=c;
        yr_vecdd(1,1)=0;
    end % c loop

    %set population size year step vector
    pop_vecdd=ones(tlimit+1,1);
    pop_vecdd(1,1)=sum(ndd);

    % set population extinction indicator vector
    ext_vecdd=zeros(tlimit+1,1);
    ext_vecdd(1,1)=0;

    %then iterate
    for iy=1:tlimit;

        if (disdd_new(iy)==1)
            ndd=adisdd*ndd;
        else
            ndd=amax*ndd;
        end

        % has the population gone extinct (i.e., n vector sums to < 1)?
        ext_thresh = 50; % set extinction threshold
        
        if sum(ndd) < ext_thresh;
            ext_vecdd(iy+1,1)=1;
        else
            ext_vecdd(iy+1,1)=0;
        end
    
        % add to population vector
        pop_vecdd(iy+1,1)=(sum(ndd));

        % Set negative density feedback function for the 'normal' matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))

        % parameters with 70000 mid-point
        acoeffn=0.1968;
        bcoeffn=2.9838;
        x0n=70000;
        y0n=0.6234;

        % Redefine survival probabilities
        s2_n = y0n+(acoeffn/(1+(sum(ndd)/x0n)^bcoeffn));
        s3_n = s2_n;
        s4_n = s2_n;
        s5_n = s2_n;
    
        % the new matrix
        amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
              s1 0 0 0 0 0
              0 s2_n 0 0 0 0
              0 0 s3_n 0 0 0
              0 0 0 s4_n 0 0
              0 0 0 0 s5_n 0];

        % Set negative density feedback function for the disease matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))
        acoeff=0.5629868;
        bcoeff=3.00328;
        x0=35000;
        y0=0.037572;

        % Redefine survival probabilities
        sc1_dd = s1;
        sc2_dd = y0+(acoeff/(1+(sum(ndd)/x0)^bcoeff));
        sc3_dd = sc2_dd;
        sc4_dd = sc2_dd;
        sc5_dd = sc2_dd;
    
        % the new disease matrix
        adisdd=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
        sc1_dd 0 0 0 0 0
        0 sc2_dd 0 0 0 0
        0 0 sc3_dd 0 0 0
        0 0 0 sc4_dd 0 0
        0 0 0 0 sc5_dd 0];

    end % for iy

    % log-transform population vector
    log_pop_vecdd=log10(pop_vecdd);

    % place population size vectors into storage matrix
    pop_vecdd_ci(:,idd) = pop_vecdd;
    log_vecdd_ci(:,idd) = log_pop_vecdd;

    % did the population go extinct during this run?
    extinct(idd,1) = sum(ext_vecdd);

end % idd loop

% Calculate mean population sizes
mean_popdd = mean(pop_vecdd_ci,2);
%mean_log_popdd = mean(log_vecdd_ci,2);
mean_log_popdd = log10(mean_popdd);

% calculate lowest percentage of original population size
pcpop_vecdd_ci = (pop_vecdd_ci/N);
minpcdd = (min(pcpop_vecdd_ci,[],1))*100;
mean_minpcdd = (mean(minpcdd));

% Calculate probability of extinction
ext = zeros(iter,1);     
for ie = 1:iter;
        if extinct(ie,1) > 0;
            ext(ie,1) = 1;
        else
            ext(ie,1) = 0;
        end
    end

pr_ext = (sum(ext)) / iter

% Calculate population confidence intervals
ord_popdd = sort(pop_vecdd_ci,2);
ord_log_popdd = sort(log_vecdd_ci,2);

% Calculate subscripts to estimate confidence intervals
fivepc = 0.05*iter;
sublo = round(fivepc/2);
subup = round(iter - sublo);

% Calculate confidence intervals for population sizes
lo_popdd = ord_popdd(:,sublo);
up_popdd = ord_popdd(:,subup);

lo_log_popdd = ord_log_popdd(:,sublo);
up_log_popdd = ord_log_popdd(:,subup);

% Plot mean trajectories & 95% confidence intervals
subplot(1,2,1), plot(yr_vecdd,mean_popdd);
axis square;
axis([-0.5 tlimit+1 0 (max(up_popdd)+(0.25*(max(mean_popdd))))]);
xlabel('year');
ylabel('N females');
hold
subplot(1,2,1), plot(yr_vecdd,up_popdd,'r:');
subplot(1,2,1), plot(yr_vecdd,lo_popdd,'r:');
hold
subplot(1,2,2), plot(yr_vecdd,mean_log_popdd);
axis square;
axis([-0.5 tlimit+1 0 (max(up_log_popdd)+(0.25*(max(mean_log_popdd))))]);
xlabel('year');
ylabel('log N females');
hold
subplot(1,2,2), plot(yr_vecdd,up_log_popdd,'r:');
subplot(1,2,2), plot(yr_vecdd,lo_log_popdd,'r:');
hold

% What is mean minimum per cent of initial population size?
mean_minpcdd

% Confidence interval
sort_minpcdd = sort(minpcdd);
minpcddsub_lo = round(0.025 * iter);
minpcddsub_hi = round(0.975 * iter);
minpcdd_lo = sort_minpcdd(minpcddsub_lo)
minpcdd_hi = sort_minpcdd(minpcddsub_hi)

subplot(1,1,1), hist(minpcdd);
axis square;
axis([0 100 0 (max(hist(minpcdd))+10)]);
xlabel('mean min % N');
ylabel('frequency');
hold


% --------------------------------------------------------------------------------------------------------------------------
% stochastic model to build disease bouts (Markhov Chain dependency) with negative density feedback (surv) - one random set
% --------------------------------------------------------------------------------------------------------------------------

% reset amax and adisease to original values
amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
       s1 0 0 0 0 0
       0 s2 0 0 0 0
       0 0 s3 0 0 0
       0 0 0 s4 0 0
       0 0 0 0 s5 0];

adisease=[0 s0*m_imm*x s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
          s1 0 0 0 0 0
          0 sc2 0 0 0 0
          0 0 sc3 0 0 0
          0 0 0 sc4 0 0
          0 0 0 0 sc5 0];

adisddd = adisease;

% initialise p_dep
p_dep = 0.50;

% set new time limit
tlimit = 200;

    % generate random number tlimit vector
    random=rand(1,tlimit);

    % sample p_disease from possible range (defined under previous model)
    p_diseaseddd = (rand * p_diseasedd_dif) + min_p_disease;

    % initial stochastic determination (1=disease; 2=no disease);
    yddd=(p_diseaseddd<=random);
    disddd=yddd+1;

    % start Markhov Chain
    disddd_new=disddd;

    for ip=1:tlimit-1;
    
        if (disddd_new(ip)==1) & (rand(1,1) <= p_dep)
             disddd_new(ip+1) = 1;
        else
             disddd_new(ip+1) = disddd_new(ip+1);
        end

    end % for ip

    % pick initial vectors
    nddd=ones(k,1);
    nddd=ssd_inc*N;
    
    %set year step vector
    yr_vecddd=ones(tlimit+1,1);
    for c=1:tlimit
        yr_vecddd(c+1,1)=c;
        yr_vecddd(1,1)=0;
    end % c loop

    %set population size year step vector
    pop_vecddd=ones(tlimit+1,1);
    pop_vecddd(1,1)=sum(nddd);

    % pick initial vector, summing to 1
    nstochddd=ones(k,1)/k;

    %then iterate
    for iy=1:tlimit;

        if (disddd_new(iy)==1)
            nddd=adisddd*nddd;
        else
            nddd=amax*nddd;
        end
    
        % add to population vector
        pop_vecddd(iy+1,1)=(sum(nddd));
    
        % one-step growth rates
        rddd(iy)=sum(nddd);
    
        % rescale vector
        nstochddd=nddd/sum(nddd);

        % Set negative density feedback function for the 'normal' matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))

        % parameters with 70000 mid-point
        acoeffn=0.1968;
        bcoeffn=2.9838;
        x0n=70000;
        y0n=0.6234;

        % Redefine survival probabilities
        s2_n = y0n+(acoeffn/(1+(sum(nddd)/x0n)^bcoeffn));
        s3_n = s2_n;
        s4_n = s2_n;
        s5_n = s2_n;
    
        % the new matrix
        amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
              s1 0 0 0 0 0
              0 s2_n 0 0 0 0
              0 0 s3_n 0 0 0
              0 0 0 s4_n 0 0
              0 0 0 0 s5_n 0];

        % Set negative density feedback function for the disease matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))
        acoeff=0.5629868;
        bcoeff=3.00328;
        x0=35000;
        y0=0.037572;

        % Redefine survival probabilities
        sc1_ddd = s1;
        sc2_ddd = y0+(acoeff/(1+(sum(nddd)/x0)^bcoeff));
        sc3_ddd = sc2_ddd;
        sc4_ddd = sc2_ddd;
        sc5_ddd = sc2_ddd;
    
        % the new disease matrix
        adisddd=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
        sc1_ddd 0 0 0 0 0
        0 sc2_ddd 0 0 0 0
        0 0 sc3_ddd 0 0 0
        0 0 0 sc4_ddd 0 0
        0 0 0 0 sc5_ddd 0];

    end % for iy

    % log-transform population vector
    log_pop_vecddd=log10(pop_vecddd);

% running mean smoother
rml=10; % number of years over which to calculate smoother
rml2=rml-1;
rmvec=zeros(tlimit-rml2,1);
		for n=1:tlimit-rml2;

            rmln=n+rml2;
			for m=n:rmln;
				rmm=mean(pop_vecddd(n:rmln,:));
            end
            
			rmvec(n,1)=rmm;
        end
 
% find subscripts to modify running mean plot
rm_sub1 = round(rml/2);

if mod(rml,2) == 1
    rm_sub2 = rm_sub1;
else
    rm_sub2 = rm_sub1+1;
end

    % Plot trajectory
    subplot(1,2,1), plot(yr_vecddd,pop_vecddd);
    axis square;
    axis([-0.5 tlimit+1 0 (max(pop_vecddd)+(0.25*(max(pop_vecddd))))]);
    xlabel('year');
    ylabel('N females');
    subplot(1,2,2), plot(yr_vecddd(rm_sub1:length(yr_vecddd)-rm_sub2,:),rmvec,'r:');
    axis square;
    axis([-0.5 tlimit+1 0 (max(pop_vecddd)+(0.25*(max(pop_vecddd))))]);
    xlabel('year');
    ylabel('N females (running mean)');

% ---------------------------------------
% plot paper figure trio (dd_surv, di_ac)
% ---------------------------------------
  
% Plot mean trajectories & 95% confidence intervals
subplot(1,3,1), plot(yr_vecdd,mean_popdd);
axis square;
axis([-0.5 tlimit+1 0 (max(up_popdd)+(0.25*(max(mean_popdd))))]);
xlabel('year');
ylabel('N females');
hold
subplot(1,3,1), plot(yr_vecdd,up_popdd,'r:');
subplot(1,3,1), plot(yr_vecdd,lo_popdd,'r:');
hold
subplot(1,3,2), plot(yr_vecddd,pop_vecddd);
axis square;
axis([-0.5 tlimit+1 0 (max(pop_vecddd)+(0.25*(max(pop_vecddd))))]);
xlabel('year');
ylabel('N females');
subplot(1,3,3), plot(yr_vecddd(rm_sub1:length(yr_vecddd)-rm_sub2,:),rmvec,'r:');
axis square;
axis([-0.5 tlimit+1 0 (max(pop_vecddd)+(0.25*(max(pop_vecddd))))]);
xlabel('year');
ylabel('N females (running mean)');

% ---------------------------------------------------------------------------------------------------------------------------------------
% stochastic model to build disease bouts (Markhov Chain dependency) with negative density feedback in survival and inter-year dependency
% ---------------------------------------------------------------------------------------------------------------------------------------

% reset amax and adisease to original values
amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
       s1 0 0 0 0 0
       0 s2 0 0 0 0
       0 0 s3 0 0 0
       0 0 0 s4 0 0
       0 0 0 0 s5 0];

adisease=[0 s0*m_imm*x s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
          s1 0 0 0 0 0
          0 sc2 0 0 0 0
          0 0 sc3 0 0 0
          0 0 0 sc4 0 0
          0 0 0 0 sc5 0];

adisddd = adisease;
      
% initialise p_dep at median value
p_dep = 0.50;

% Re-define tlimit
tlimit = 200;

% define matrix to hold population sizes
pop_vecddd_ci = ones(tlimit+1,iter);
log_vecddd_ci = ones(tlimit+1,iter);

% set extinction sum vector
  extinct=zeros(iter,1);

% do loop to estimate population size confidence intervals
for iddd=1:iter;

    % generate random number tlimit vector
    random=rand(1,tlimit);

    % sample p_disease from possible range (defined under previous model)
    p_diseaseddd = (rand * p_diseasedd_dif) + min_p_disease;

    % initial stochastic determination (1=disease; 2=no disease);
    yddd=(p_diseaseddd<=random);
    disddd=yddd+1;

    % start Markhov Chain
    disddd_new=disddd;

    for ip=1:tlimit-1;
    
        if (disddd_new(ip)==1) & (rand(1,1) <= p_dep)
             disddd_new(ip+1) = 1;
        else
             disddd_new(ip+1) = disddd_new(ip+1);
        end

    end % for ip

    % pick initial vectors
    nddd=ones(k,1);
    nddd=ssd_inc*N;
    
    %set year step vector
    yr_vecddd=ones(tlimit+1,1);
    for c=1:tlimit
        yr_vecddd(c+1,1)=c;
        yr_vecddd(1,1)=0;
    end % c loop

    % set population extinction indicator vector
    ext_vecm=zeros(tlimit+1,1);
    ext_vecm(1,1)=0;

    %set population size year step vector
    pop_vecddd=ones(tlimit+1,1);
    pop_vecddd(1,1)=sum(nddd);

    % pick initial vector, summing to 1
    nstochddd=ones(k,1)/k;

    %then iterate
    for iy=1:tlimit;

        if (disddd_new(iy)==1)
            nddd=adisddd*nddd;
        else
            nddd=amax*nddd;
        end

        % has the population gone extinct (i.e., n vector sums to < 1)?
        ext_thresh = 50; % set extinction threshold
        
        if sum(nddd) < ext_thresh;
            ext_vecm(iy+1,1)=1;
        else
            ext_vecm(iy+1,1)=0;
        end
    
        % add to population vector
        pop_vecddd(iy+1,1)=(sum(nddd));
    
        % one-step growth rates
        rddd(iy)=sum(nddd);
    
        % rescale vector
        nstochddd=nddd/sum(nddd);

        % Set negative density feedback function for the 'normal' matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))

        % parameters with 70000 mid-point
        acoeffn=0.1968;
        bcoeffn=2.9838;
        x0n=70000;
        y0n=0.6234;

        % Redefine survival probabilities
        s2_n = y0n+(acoeffn/(1+(sum(nddd)/x0n)^bcoeffn));
        s3_n = s2_n;
        s4_n = s2_n;
        s5_n = s2_n;
    
        % the new matrix
        amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
              s1 0 0 0 0 0
              0 s2_n 0 0 0 0
              0 0 s3_n 0 0 0
              0 0 0 s4_n 0 0
              0 0 0 0 s5_n 0];

        % Set negative density feedback function for the disease matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))
        acoeff=0.5629868;
        bcoeff=3.00328;
        x0=35000;
        y0=0.037572;

        % Redefine survival probabilities
        sc1_ddd = s1;
        sc2_ddd = y0+(acoeff/(1+(sum(nddd)/x0)^bcoeff));
        sc3_ddd = sc2_ddd;
        sc4_ddd = sc2_ddd;
        sc5_ddd = sc2_ddd;
    
        % the new disease matrix
        adisddd=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
        sc1_ddd 0 0 0 0 0
        0 sc2_ddd 0 0 0 0
        0 0 sc3_ddd 0 0 0
        0 0 0 sc4_ddd 0 0
        0 0 0 0 sc5_ddd 0];

        % define negative density feedback on probability of inter-year disease dependency
        % Re-define 4-parameter logistic function parameters
        acoeffdep=-1.0944513806541;
        bcoeffdep=3.39495801978165;
        x0dep=35000;
        y0dep=1.09501076129535;
        
        p_dep = y0dep+(acoeffdep/(1+(sum(ndd)/x0dep)^bcoeffdep));
        
    end % for iy

    % log-transform population vector
    log_pop_vecddd=log10(pop_vecddd);

    % place population size vectors into storage matrix
    pop_vecddd_ci(:,iddd) = pop_vecddd;
    log_vecddd_ci(:,iddd) = log_pop_vecddd;

end % iddd loop

% Calculate mean population sizes
mean_popddd = mean(pop_vecddd_ci,2);
%mean_log_popddd = mean(log_vecddd_ci,2);
mean_log_popddd = log10(mean_popddd);

% calculate lowest percentage of original population size
pcpop_vecddd_ci = (pop_vecddd_ci/N);
minpcddd = (min(pcpop_vecddd_ci,[],1))*100;
mean_minpcddd = (mean(minpcddd));

% Calculate probability of extinction
ext = zeros(iter,1);     
for ie = 1:iter;
        if extinct(ie,1) > 0;
            ext(ie,1) = 1;
        else
            ext(ie,1) = 0;
        end
    end

pr_ext = (sum(ext)) / iter

% Calculate population confidence intervals
ord_popddd = sort(pop_vecddd_ci,2);
ord_log_popddd = sort(log_vecddd_ci,2);

% Calculate subscripts to estimate confidence intervals
fivepc = 0.05*iter;
sublo = round(fivepc/2);
subup = round(iter - sublo);

% Calculate confidence intervals for population sizes
lo_popddd = ord_popddd(:,sublo);
up_popddd = ord_popddd(:,subup);

lo_log_popddd = ord_log_popddd(:,sublo);
up_log_popddd = ord_log_popddd(:,subup);

% Plot mean trajectories & 95% confidence intervals
subplot(1,2,1), plot(yr_vecddd,mean_popddd);
axis square;
axis([-0.5 tlimit+1 0 (max(up_popddd)+(0.25*(max(mean_popddd))))]);
xlabel('year');
ylabel('N females');
hold
subplot(1,2,1), plot(yr_vecddd,up_popddd,'r:');
subplot(1,2,1), plot(yr_vecddd,lo_popddd,'r:');
hold
subplot(1,2,2), plot(yr_vecddd,mean_log_popddd);
axis square;
axis([-0.5 tlimit+1 0 (max(up_log_popddd)+(0.25*(max(mean_log_popddd))))]);
xlabel('year');
ylabel('log N females');
hold
subplot(1,2,2), plot(yr_vecddd,up_log_popddd,'r:');
subplot(1,2,2), plot(yr_vecddd,lo_log_popddd,'r:');
hold

% What is mean minimum per cent of initial population size?
mean_minpcddd

% Estimate CI
sort_minpcddd = sort(minpcddd);
minpcdddsub_lo = round(0.025 * iter);
minpcdddsub_hi = round(0.975 * iter);
minpcddd_lo = sort_minpcdd(minpcdddsub_lo)
minpcddd_hi = sort_minpcdd(minpcdddsub_hi)

subplot(1,1,1), hist(minpcddd);
axis square;
axis([0 100 0 (max(hist(minpcddd))+10)]);
xlabel('mean min % N');
ylabel('frequency');
hold

% ---------------------------------------------------------------------------------------------------------------------------------
% stochastic model to build disease bouts (Markhov Chain dependency) with negative density feedback (surv & p_dep) - one random set
% ---------------------------------------------------------------------------------------------------------------------------------

% reset amax and adisease to original values
amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
       s1 0 0 0 0 0
       0 s2 0 0 0 0
       0 0 s3 0 0 0
       0 0 0 s4 0 0
       0 0 0 0 s5 0];

adisease=[0 s0*m_imm*x s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
          s1 0 0 0 0 0
          0 sc2 0 0 0 0
          0 0 sc3 0 0 0
          0 0 0 sc4 0 0
          0 0 0 0 sc5 0];

adisddd = adisease;

% initialise p_dep
p_dep = 0.50;

% set new time limit
tlimit = 200;

    % generate random number tlimit vector
    random=rand(1,tlimit);

    % sample p_disease from possible range (defined under previous model)
    p_diseaseddd = (rand * p_diseasedd_dif) + min_p_disease;

    % initial stochastic determination (1=disease; 2=no disease);
    yddd=(p_diseaseddd<=random);
    disddd=yddd+1;

    % start Markhov Chain
    disddd_new=disddd;

    for ip=1:tlimit-1;
    
        if (disddd_new(ip)==1) & (rand(1,1) <= p_dep)
             disddd_new(ip+1) = 1;
        else
             disddd_new(ip+1) = disddd_new(ip+1);
        end

    end % for ip

    % pick initial vectors
    nddd=ones(k,1);
    nddd=ssd_inc*N;
    
    %set year step vector
    yr_vecddd=ones(tlimit+1,1);
    for c=1:tlimit
        yr_vecddd(c+1,1)=c;
        yr_vecddd(1,1)=0;
    end % c loop

    %set population size year step vector
    pop_vecddd=ones(tlimit+1,1);
    pop_vecddd(1,1)=sum(nddd);

    % pick initial vector, summing to 1
    nstochddd=ones(k,1)/k;

    %then iterate
    for iy=1:tlimit;

        if (disddd_new(iy)==1)
            nddd=adisddd*nddd;
        else
            nddd=amax*nddd;
        end
    
        % add to population vector
        pop_vecddd(iy+1,1)=(sum(nddd));
    
        % one-step growth rates
        rddd(iy)=sum(nddd);
    
        % rescale vector
        nstochddd=nddd/sum(nddd);

        % Set negative density feedback function for the 'normal' matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))

        % parameters with 70000 mid-point
        acoeffn=0.1968;
        bcoeffn=2.9838;
        x0n=70000;
        y0n=0.6234;

        % Redefine survival probabilities
        s2_n = y0n+(acoeffn/(1+(sum(nddd)/x0n)^bcoeffn));
        s3_n = s2_n;
        s4_n = s2_n;
        s5_n = s2_n;
    
        % the new matrix
        amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
              s1 0 0 0 0 0
              0 s2_n 0 0 0 0
              0 0 s3_n 0 0 0
              0 0 0 s4_n 0 0
              0 0 0 0 s5_n 0];

        % Set negative density feedback function for the disease matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))
        acoeff=0.5629868;
        bcoeff=3.00328;
        x0=35000;
        y0=0.037572;

        % Redefine survival probabilities
        sc1_ddd = s1;
        sc2_ddd = y0+(acoeff/(1+(sum(nddd)/x0)^bcoeff));
        sc3_ddd = sc2_ddd;
        sc4_ddd = sc2_ddd;
        sc5_ddd = sc2_ddd;
    
        % the new disease matrix
        adisddd=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
        sc1_ddd 0 0 0 0 0
        0 sc2_ddd 0 0 0 0
        0 0 sc3_ddd 0 0 0
        0 0 0 sc4_ddd 0 0
        0 0 0 0 sc5_ddd 0];

        % define negative density feedback on probability of inter-year disease dependency
        % Re-define 4-parameter logistic function parameters
        acoeffdep=-0.8771;
        bcoeffdep=3.3825;
        x0dep=35000;
        y0dep=0.9767;
        
        p_dep = y0dep+(acoeffdep/(1+(sum(nddd)/x0dep)^bcoeffdep));

    end % for iy

    % log-transform population vector
    log_pop_vecddd=log10(pop_vecddd);

% find periodicity of trend

% running mean smoother
rml=10; % number of years over which to calculate smoother
rml2=rml-1;
rmvec=zeros(tlimit-rml2,1);
		for n=1:tlimit-rml2;

            rmln=n+rml2;
			for m=n:rmln;
				rmm=mean(pop_vecddd(n:rmln,:));
            end
            
			rmvec(n,1)=rmm;
        end
 
% find subscripts to modify running mean plot
rm_sub1 = round(rml/2);

if mod(rml,2) == 1
    rm_sub2 = rm_sub1;
else
    rm_sub2 = rm_sub1+1;
end

    % Plot trajectory
    subplot(1,2,1), plot(yr_vecddd,pop_vecddd);
    axis square;
    axis([-0.5 tlimit+1 0 (max(pop_vecddd)+(0.25*(max(pop_vecddd))))]);
    xlabel('year');
    ylabel('N females');
    subplot(1,2,2), plot(yr_vecddd(rm_sub1:length(yr_vecddd)-rm_sub2,:),rmvec,'r:');
    axis square;
    axis([-0.5 tlimit+1 0 (max(pop_vecddd)+(0.25*(max(pop_vecddd))))]);
    xlabel('year');
    ylabel('N females (running mean)');
        
       
% calculate differences between successive population sizes
pop_dif = diff(rmvec);
pop_dir = sign(pop_dif);
pop_difr = diff(pop_dir);

% set counter vectors
neg_tr_num=0;
pos_tr_num=0;

% calculate number of positive and negative trends in time series
    for d=1:(length(pop_dir));
        
        if pop_dir(d,:) < 0;
            neg_tr_num = neg_tr_num+1;
        else
            pos_tr_num = pos_tr_num+1;
        end
        
    end

% Count bouts of increase & decline
neg_sub = find(pop_difr==-2);
pos_sub = find(pop_difr==2);
start_trend = pop_dir(1,:);
lneg = length(neg_sub);
lpos = length(pos_sub);
bouts = zeros(lneg+lpos,1);
bouts(1:lneg,:)=neg_sub;
bouts(lneg+1:lneg+lpos,:)=pos_sub;
bouts=sort(bouts);
bouts_dif = diff(bouts);
bouts_diff = zeros(length(bouts)+1,1);
bouts_diff(1,:)=bouts(1,:);
bouts_diff(2:length(bouts),1)=bouts_dif;
bouts_diff(length(bouts_diff),1)=length(pop_dir) - bouts(length(bouts)) + 1;

% Get neg and pos bouts in separate vectors
neg_bouts=zeros(lneg,1);
pos_bouts=zeros(lpos,1);

if (length(bouts_diff)/2) - round(length(bouts_diff)/2) > 0;
    addd = 1;
else
    addd = 0;
end

length_neg = lneg;
length_pos = lpos;

if start_trend == -1;
    length_neg = length_neg + 1;
else
    length_pos = length_pos + 1;
end


if start_trend == -1;
    for n = 1:2:length(bouts_diff);
        neg_bouts(n,1)=bouts_diff(n,:);
    end
    for p = 2:2:length(bouts_diff);
        pos_bouts(p,1)=bouts_diff(p,:);
    end
    negbouts = neg_bouts(find(neg_bouts>0),1);
    posbouts = pos_bouts(find(pos_bouts>0),1);
else
    for n = 1:2:length(bouts_diff);
        pos_bouts(n,1)=bouts_diff(n,:);
    end
    for p = 2:2:length(bouts_diff);
        neg_bouts(p,1)=bouts_diff(p,:);
    end
    negbouts = neg_bouts(find(neg_bouts>0),1);
    posbouts = pos_bouts(find(pos_bouts>0),1);
end

% calculate mean (& SD) bout lengths (years)
negbout_m = mean(negbouts)
negbout_sd = std(negbouts)

posbout_m = mean(posbouts)
posbout_sd = std(posbouts)

% ---------------------------------------
% plot paper figure trio (dd_surv, dd_ac)
% ---------------------------------------
  
% Plot mean trajectories & 95% confidence intervals
subplot(1,3,1), plot(yr_vecddd,mean_popddd);
axis square;
axis([-0.5 tlimit+1 0 (max(up_popddd)+(0.25*(max(mean_popddd))))]);
xlabel('year');
ylabel('N females');
hold
subplot(1,3,1), plot(yr_vecddd,up_popddd,'r:');
subplot(1,3,1), plot(yr_vecddd,lo_popddd,'r:');
hold
subplot(1,3,2), plot(yr_vecddd,pop_vecddd);
axis square;
axis([-0.5 tlimit+1 0 (max(pop_vecddd)+(0.25*(max(pop_vecddd))))]);
xlabel('year');
ylabel('N females');
subplot(1,3,3), plot(yr_vecddd(rm_sub1:length(yr_vecddd)-rm_sub2,:),rmvec,'r:');
axis square;
axis([-0.5 tlimit+1 0 (max(pop_vecddd)+(0.25*(max(pop_vecddd))))]);
xlabel('year');
ylabel('N females (running mean)');

% ---------------------------------------------------------------------------------------------------------------------------------------
% stochastic model to build disease bouts (Markhov Chain dependency) with negative density feedback (surv & p_dep) - estimate periodicity
% ---------------------------------------------------------------------------------------------------------------------------------------

% reset amax and adisease to original values
amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
       s1 0 0 0 0 0
       0 s2 0 0 0 0
       0 0 s3 0 0 0
       0 0 0 s4 0 0
       0 0 0 0 s5 0];

adisease=[0 s0*m_imm*x s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
          s1 0 0 0 0 0
          0 sc2 0 0 0 0
          0 0 sc3 0 0 0
          0 0 0 sc4 0 0
          0 0 0 0 sc5 0];

adisddd = adisease;

% initialise p_dep
p_dep = 0.50;

% set new time limit
tlimit = 200;

% run i iterations to estimate average (& variance) duration (years) of positive population
% trends and negative population trends

iter = 1000; % set number of iterations

negbout_mm = zeros(iter,1);
posbout_mm = zeros(iter,1);

hertz1_vec = zeros(iter,1);
hertz2_vec = zeros(iter,1);

for i = 1:iter;

    % generate random number tlimit vector
    random=rand(1,tlimit);

    % sample p_disease from possible range (defined under previous model)
    p_diseaseddd = (rand * p_diseasedd_dif) + min_p_disease;

    % initial stochastic determination (1=disease; 2=no disease);
    yddd=(p_diseaseddd<=random);
    disddd=yddd+1;

    % start Markhov Chain
    disddd_new=disddd;

    for ip=1:tlimit-1;
    
        if (disddd_new(ip)==1) & (rand(1,1) <= p_dep)
             disddd_new(ip+1) = 1;
        else
             disddd_new(ip+1) = disddd_new(ip+1);
        end

    end % for ip

    % pick initial vectors
    nddd=ones(k,1);
    %nddd(1:2,1)=0.2382*N/2; % Guiler 1978
    %nddd(3:6,1)=0.7618*N/4; % Guiler 1978
    nddd=ssd_inc*N;
    
    %set year step vector
    yr_vecddd=ones(tlimit+1,1);
    for c=1:tlimit
        yr_vecddd(c+1,1)=c;
        yr_vecddd(1,1)=0;
    end % c loop

    %set population size year step vector
    pop_vecddd=ones(tlimit+1,1);
    pop_vecddd(1,1)=sum(nddd);

    % pick initial vector, summing to 1
    nstochddd=ones(k,1)/k;

    %then iterate
    for iy=1:tlimit;

        if (disddd_new(iy)==1)
            nddd=adisddd*nddd;
        else
            nddd=amax*nddd;
        end
    
        % add to population vector
        pop_vecddd(iy+1,1)=(sum(nddd));
    
        % one-step growth rates
        rddd(iy)=sum(nddd);
    
        % rescale vector
        nstochddd=nddd/sum(nddd);

        % Set negative density feedback function for the 'normal' matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))

        % parameters with 70000 mid-point
        acoeffn=0.1968;
        bcoeffn=2.9838;
        x0n=70000;
        y0n=0.6234;

        % Redefine survival probabilities
        s2_n = y0n+(acoeffn/(1+(sum(nddd)/x0n)^bcoeffn));
        s3_n = s2_n;
        s4_n = s2_n;
        s5_n = s2_n;
    
        % the new matrix
        amax=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
              s1 0 0 0 0 0
              0 s2_n 0 0 0 0
              0 0 s3_n 0 0 0
              0 0 0 s4_n 0 0
              0 0 0 0 s5_n 0];

        % Set negative density feedback function for the disease matrix (f(x)=y0+(acoeff/(1+(x/x0)^bcoeff)))
        acoeff=0.5629868;
        bcoeff=3.00328;
        x0=35000;
        y0=0.037572;

        % Redefine survival probabilities
        sc1=s1;
        sc2_ddd = y0+(acoeff/(1+(sum(nddd)/x0)^bcoeff));
        sc3_ddd = sc2_ddd;
        sc4_ddd = sc2_ddd;
        sc5_ddd = sc2_ddd;
        
        % the new disease matrix
        adisddd=[s0*m_imm*x s0*m_primi*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_prime*x*rs s0*m_old*x*rs
        sc1_ddd 0 0 0 0 0
        0 sc2_ddd 0 0 0 0
        0 0 sc3_ddd 0 0 0
        0 0 0 sc4_ddd 0 0
        0 0 0 0 sc5_ddd 0];

        % define negative density feedback on probability of inter-year disease dependency
        % Re-define 4-parameter logistic function parameters
        acoeffdep=-0.8771;
        bcoeffdep=3.3825;
        x0dep=35000;
        y0dep=0.9767;
        
        p_dep = y0dep+(acoeffdep/(1+(sum(nddd)/x0dep)^bcoeffdep));

    end % for iy

    % log-transform population vector
    log_pop_vecddd=log10(pop_vecddd);

    % find periodicity of trend
    % running mean smoother
    rml=10; % number of years over which to calculate smoother
    rml2=rml-1;
    rmvec=zeros(tlimit-rml2,1);
		    for n=1:tlimit-rml2;

              rmln=n+rml2;
    			for m=n:rmln;
				    rmm=mean(pop_vecddd(n:rmln,:));
                end
            
			    rmvec(n,1)=rmm;
           end

        % find subscripts to modify running mean plot
    rm_sub1 = round(rml/2);

    if mod(rml,2) == 1
        rm_sub2 = rm_sub1;
    else
        rm_sub2 = rm_sub1+1;
    end

    % calculate differences between successive population sizes
    pop_dif = diff(rmvec);
    pop_dir = sign(pop_dif);
    pop_difr = diff(pop_dir);

    % set counter vectors
    neg_tr_num=0;
    pos_tr_num=0;

    % calculate number of positive and negative trends in time series
    for d=1:(length(pop_dir));
        
        if pop_dir(d,:) < 0;
            neg_tr_num = neg_tr_num+1;
        else
            pos_tr_num = pos_tr_num+1;
        end
        
    end

    % Count bouts of increase & decline
    neg_sub = find(pop_difr==-2);
    pos_sub = find(pop_difr==2);
    start_trend = pop_dir(1,:);
    lneg = length(neg_sub);
    lpos = length(pos_sub);
    bouts = zeros(lneg+lpos,1);
    bouts(1:lneg,:)=neg_sub;
    bouts(lneg+1:lneg+lpos,:)=pos_sub;
    bouts=sort(bouts);
    bouts_dif = diff(bouts);
    bouts_diff = zeros(length(bouts)+1,1);
    bouts_diff(1,:)=bouts(1,:);
    bouts_diff(2:length(bouts),1)=bouts_dif;
    bouts_diff(length(bouts_diff),1)=length(pop_dir) - bouts(length(bouts)) + 1;

    % Get neg and pos bouts in separate vectors
    neg_bouts=zeros(lneg,1);
    pos_bouts=zeros(lpos,1);

    if (length(bouts_diff)/2) - round(length(bouts_diff)/2) > 0;
        addd = 1;
    else
        addd = 0;
    end

    length_neg = lneg;
    length_pos = lpos;

    if start_trend == -1;
        length_neg = length_neg + 1;
    else
        length_pos = length_pos + 1;
    end

    if start_trend == -1;
        for n = 1:2:length(bouts_diff);
            neg_bouts(n,1)=bouts_diff(n,:);
        end
        for p = 2:2:length(bouts_diff);
            pos_bouts(p,1)=bouts_diff(p,:);
        end
        negbouts = neg_bouts(find(neg_bouts>0),1);
        posbouts = pos_bouts(find(pos_bouts>0),1);
    else
        for n = 1:2:length(bouts_diff);
            pos_bouts(n,1)=bouts_diff(n,:);
        end
        for p = 2:2:length(bouts_diff);
            neg_bouts(p,1)=bouts_diff(p,:);
        end
        negbouts = neg_bouts(find(neg_bouts>0),1);
        posbouts = pos_bouts(find(pos_bouts>0),1);
    end

    % calculate mean bout lengths (years)
    negbout_m = mean(negbouts);
    posbout_m = mean(posbouts);

    % place means over time interval in iteration mean vector
    negbout_mm(i,1) = negbout_m;
    posbout_mm(i,1) = posbout_m;

    i

    % fast Fourier transform (FTT) analysis to estimate frequency of major bouts within time interval
    spvec = pop_vecddd/(max(pop_vecddd));
    Yf = fft(spvec);
    Nf = length(spvec);
    Tf = 1;
    tf = [0:Nf-1]/Nf;
    tf = tf*Tf;
    pf = abs(fft(spvec))/(Nf/2);
    pf = pf(1:Nf/2).^2;
    freq = [0:Nf/2-1]/Tf;
    
    % first major cycle
    hertz1 = (find(pf==max(pf(2:10)))) - 1;
    ep_int1 = tlimit/hertz1;

    % second major cycle
    hertz2 = (find(pf==max(pf(3:10)))) - 1;
    ep_int2 = tlimit/hertz2;

    % Place cycle durations into vectors
    epint1_vec(i,1) = ep_int1;
    epint2_vec(i,1) = ep_int2;
    
end % i loop    

% calculate means and CI of major cycle durations over all iterations
epint1_lo = quantile(epint1_vec,0.025)
epint1_mean = mean(epint1_vec)
epint1_up = quantile(epint1_vec,0.975)

epint2_lo = quantile(epint2_vec,0.025)
epint2_mean = mean(epint2_vec)
epint2_up = quantile(epint2_vec,0.975)

% plot histograms of major cycle durations
subplot(1,2,1), hist(epint1_vec);
axis square;
xlabel('main cycle duration (years)');
ylabel('frequency')
subplot(1,2,2), hist(epint2_vec);
axis square;
xlabel('secondary cycle duration (years)');
ylabel('frequency')

pause

% calculate means and CI of periodicities of iterations
negbout_mean = mean(negbout_mm);
posbout_mean = mean(posbout_mm);

negbout_up = quantile(negbout_mm,0.975);
negbout_lo = quantile(negbout_mm,0.025);

posbout_up = quantile(posbout_mm,0.975);
posbout_lo = quantile(posbout_mm,0.025);

negbout_lod = negbout_mean - negbout_lo;
negbout_upd = negbout_up - negbout_mean;

posbout_lod = posbout_mean - posbout_lo;
posbout_upd = posbout_up - posbout_mean;

% plot average & CIs of periodicities
bout_means = zeros(2,1);
bout_means(1,:) = negbout_mean;
bout_means(2,:) = posbout_mean;

bouts_lo = zeros(2,1);
bouts_lo(1,:) = negbout_lod;
bouts_lo(2,:) = posbout_lod;

bouts_hi = zeros(2,1);
bouts_up(1,:) = negbout_upd;
bouts_up(2,:) = posbout_upd;

% Plot error bar graph
errorbar(1:2,bout_means,bouts_lo,bouts_up);
hold
bar(1:2,bout_means,0.3);
xlabel('trend (1=neg; 2=pos)');
ylabel('duration (years)')
axis square;
axis([0.5 2.5 0 max(negbout_up,posbout_up)+5]);
hold

bout_means
bout_means - bouts_lo
bout_means + bouts_up
