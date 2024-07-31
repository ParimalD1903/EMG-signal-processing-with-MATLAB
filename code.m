%% Seperation of grasps and repetition and storing it in one cell
id= dir;
 for j = 3:length(id)
     k=j-2;
data = load(id(j).name);
emg_data = {};
for n = 1:23
a=data.stimulus;
if(a(a>=1))
a(a~=n)=0;
a(a==n)=1;
end
for m = 1:6
    c = a.*(data.repetition);
    if (m>=1)
        c(c~=m)=0;
        c(c==m)=1;
        e = single(c);
        d =e.*(data.emg); 
        d = double(d);
        
    end
    x = any(d==0,2);
    d(x,:) = [];
    emg_rep{m}=d;
end
emg_sub{k,n} = emg_rep;
end
disp(k)
 end


%% Selection of grasps

 a = emg_sub;
a(:,2) = [];
a(:,7) = [];
a(:,9) = [];
a(:,11) = [];
a(:,13) = [];
a(:,14) = [];
a(:,15) = [];
a(:,15) = [];
a(:,15) = [];


%% Shifting of grasps according to first 7 being power(0) and the next 7 being precision(1) 
a(:,[6,8])=a(:,[8,6]);
a(:,[7,14])=a(:,[14,7]);

%% Windowing

frame_cell = {};
for i = 1:10
for j = 1:14
    for k = 1:6
        [r,c] = size(a{i,j}{:,k});
n_start = 1;
window_size = 500;
overlap_size = 250;
while (n_start+window_size<=r)
    window = a{i,j}{:,k}(n_start:n_start+window_size-1,:);
    frame_cell{end+1} = window;
    window_emg{i,j}{:,k} = frame_cell;
    n_start = n_start + window_size - overlap_size;
end
frame_cell(:) = [];
    end
end
end

%% Decomposition 

maxAlpha = 2000;
init_omega = 0;
fs = 100;
stopc = 3;
tol = 1e-7;
tau = 0;
kdm = 4;
emg_decom = {};
for m = 6:10
    for n = 1:14
        for p = 1:6
            for o = 1:length(window_emg{m,n}{1,p})
            d = window_emg{m,n}{1,p}{1,o};
            [u{1,p}{1,o},omega,u_hat]=ssmvmd(d,maxAlpha,init_omega,fs,stopc,tol,tau,kdm);
            u_decom{m,n} = u;
            end
        end
    end
    disp(m) 
end

%% squeeze
for m = 1:5
for n = 1:14
    for p = 1:6
        for l = 1:length(u_decom_1_5{m, n}{1, p})
            for k = 1:12
            emg_sqz1{m, n}{1, p}{1, l}{:, k}=squeeze(u_decom_1_5{m, n}{1, p}{1, l}(:,k,:));
            end
        end
    end
end
disp(m)
end

%% All EMG feature's extraction
emg_feature = {};
for m = 1:5
    for n = 1:14
        for p = 1:6
             for z = 1:length(emg_sqz1{m,n}{1,p})
                 d = cell2mat(emg_sqz1{m,n}{1,p}{1,z});
                 for o = 1:48
                     a = d(:,o);
                     MAD {m,n}{1,p}{z,o}= jMeanAbsoluteDeviation(a);
                     EMAV{m,n}{1,p}{z,o} = jEnhancedMeanAbsoluteValue(a);
                     ASS{m,n}{1,p}{z,o} = jAbsoluteValueOfTheSummationOfSquareRoot(a);
                     ASM{m,n}{1,p}{z,o} = jAbsoluteValueOfTheSummationOfExpRoot(a); 
                     MMAV{m,n}{1,p}{z,o}= jModifiedMeanAbsoluteValue(a);
                     WL{m,n}{1,p}{z,o} = jWaveformLength(a);             
                     VAR{m,n}{1,p}{z,o} = jVariance(a);
                     SD{m,n}{1,p}{z,o} = jStandardDeviation(a);
                     SKEW{m,n}{1,p}{z,o} = jSkewness(a);
                     SSI{m,n}{1,p}{z,o} = jSimpleSquareIntegral(a);
                     RMS{m,n}{1,p}{z,o} = jRootMeanSquare(a);
                     FZC{m,n}{1,p}{z,o} = jNewZeroCrossing(a);
                     MMAV2{m,n}{1,p}{z,o} = jModifiedMeanAbsoluteValue2(a);
                     MSR{m,n}{1,p}{z,o}= jMeanValueOfTheSquareRoot(a);
                     MAV{m,n}{1,p}{z,o} = jMeanAbsoluteValue(a);
                     MAD {m,n}{1,p}{z,o}= jMeanAbsoluteDeviation(a);
                     MFL{m,n}{1,p}{z,o} = jMaximumFractalLength(a);
                     LTKEO{m,n}{1,p}{z,o} = jLogTeagerKaiserEnergyOperator(a);
                     LDASDV{m,n}{1,p}{z,o} = jLogDifferenceAbsoluteStandardDeviationValue(a);
                     LDAMV{m,n}{1,p}{z,o} = jLogDifferenceAbsoluteMeanValue(a);
                     LD{m,n}{1,p}{z,o} = jLogDetector(a);
                     LCOV{m,n}{1,p}{z,o} = jLogCoefficientOfVariation(a);
                     IQR{m,n}{1,p}{z,o} = jInterquartileRange(a);
                     IEMG{m,n}{1,p}{z,o} = jIntegratedEMG(a);
                     EWL{m,n}{1,p}{z,o} = jEnhancedWaveLength(a);
                     EMAV{m,n}{1,p}{z,o} = jEnhancedMeanAbsoluteValue(a);
                     DVARV{m,n}{1,p}{z,o} = jDifferenceVarianceValue(a);
                     DASDV{m,n}{1,p}{z,o} = jDifferenceAbsoluteStandardDeviationValue(a);
                     DAMV{m,n}{1,p}{z,o} = jDifferenceAbsoluteMeanValue(a);
                     COV{m,n}{1,p}{z,o} = jCoefficientOfVariation(a);
                     ME {m,n}{1,p}{z,o}= jAverageEnergy(a);
                     AAC{m,n}{1,p}{z,o} = jAverageAmplitudeChange(a);
                     KURT{m,n}{1,p}{z,o} = jKurtosis(a);
                     ASS{m,n}{1,p}{z,o} = jAbsoluteValueOfTheSummationOfSquareRoot(a);
                     ASM{m,n}{1,p}{z,o} = jAbsoluteValueOfTheSummationOfExpRoot(a); 
                 end 
             end
        end     
    end 
    disp(m)
end 
%% feature cat
for m = 1:5
    for n = 1:14
        emav{m, n} = vertcat(EMAV{m, n}{:});
        mmav{m, n} = vertcat(MMAV{m, n}{:});
        ass{m, n} = vertcat(ASS{m, n}{:});
        asm{m, n} = vertcat(ASM{m, n}{:});
        mad{m, n} = vertcat(MAD{m, n}{:});
        aac{m, n} = vertcat(AAC{m, n}{:});
        cov{m, n} = vertcat(COV{m, n}{:});
        damv{m, n} = vertcat(DAMV{m, n}{:});
        dasdv{m, n} = vertcat(DASDV{m, n}{:});
        dvarv{m, n} = vertcat(DVARV{m, n}{:});
        ewl{m, n} = vertcat(EWL{m, n}{:});
        fzc{m, n} = vertcat(FZC{m, n}{:});
        iemg{m, n} = vertcat(IEMG{m, n}{:});
        iqr{m, n} = vertcat(IQR{m, n}{:});
        kurt{m, n} = vertcat(KURT{m, n}{:});
        lcov{m, n} = vertcat(LCOV{m, n}{:});
        ld{m, n} = vertcat(LD{m, n}{:});
        ldamv{m, n} = vertcat(LDAMV{m, n}{:});
        ldasdv{m, n} = vertcat(LDASDV{m, n}{:});
        ltkeo{m, n} = vertcat(LTKEO{m, n}{:});
        mav{m, n} = vertcat(MAV{m, n}{:});
        me{m, n} = vertcat(ME{m, n}{:});
        mfl{m, n} = vertcat(MFL{m, n}{:});
        mmav21{m, n} = vertcat(MMAV2{m, n}{:});
        msr{m, n} = vertcat(MSR{m, n}{:});
        rms{m, n} = vertcat(RMS{m, n}{:});
        sd{m, n} = vertcat(SD{m, n}{:});
        skew{m, n} = vertcat(SKEW{m, n}{:});
        ssi{m, n} = vertcat(SSI{m, n}{:});
        var{m, n} = vertcat(VAR{m, n}{:});
        wl{m, n} = vertcat(WL{m, n}{:});
 
        

    end
end

 %% feature vertcat
aacvt = vertcat(cov{:});
asmvt = vertcat(asm{:});
assvt = vertcat(ass{:});
madvt = vertcat(mad{:});
emavvt = vertcat(emav{:});
mmavvt = vertcat(mmav{:});
covvt = vertcat(cov{:});
damvvt = vertcat(damv{:});
dasdvvt = vertcat(dasdv{:});
dvarvvt = vertcat(dvarv{:});
ewlvt = vertcat(ewl{:});
fzcvt = vertcat(fzc{:});
iemgvt = vertcat(iemg{:}); 
iqrvt = vertcat(iqr{:});
kurtvt = vertcat(kurt{:});
lcovvt = vertcat(lcov{:});
ldvt = vertcat(ld{:});
ldamvvt = vertcat(ldamv{:});
ldasdvvt = vertcat(ldasdv{:});
ltkeovt = vertcat(ltkeo{:});
mavvt = vertcat(mav{:});
mevt = vertcat(me{:});
mflvt = vertcat(mfl{:});
mmav21vt = vertcat(mmav21{:});
msrvt = vertcat(msr{:});
rmsvt = vertcat(rms{:});
sdvt = vertcat(sd{:});
skewvt = vertcat(skew{:});
ssivt = vertcat(ssi{:});
varvt = vertcat(var{:});
wlvt = vertcat(wl{:});


%% feature extraction of selected features and vertical concatenation
emg_feature = {};
for m = 1:5
    for n = 1:14
        for p = 1:6
             for z = 1:length(emg_sqz1{m,n}{1,p})
                 d = cell2mat(emg_sqz1{m,n}{1,p}{1,z});
                 for o = 1:48
                     a = d(:,o);
                     MAD {m,n}{1,p}{z,o}= jMeanAbsoluteDeviation(a);
                     ASS{m,n}{1,p}{z,o} = jAbsoluteValueOfTheSummationOfSquareRoot(a);
                     ASM{m,n}{1,p}{z,o} = jAbsoluteValueOfTheSummationOfExpRoot(a); 
                     SD{m,n}{1,p}{z,o} = jStandardDeviation(a);
                     IEMG{m,n}{1,p}{z,o} = jIntegratedEMG(a);
                     RMS{m,n}{1,p}{z,o} = jRootMeanSquare(a);
                     MAV{m,n}{1,p}{z,o} = jMeanAbsoluteValue(a);

                    
                 end 
             end

        end 
     ass{m, n} = vertcat(ASS{m, n}{:});
     asm{m, n} = vertcat(ASM{m, n}{:});
     mad{m, n} = vertcat(MAD{m, n}{:});
     rms{m, n} = vertcat(RMS{m, n}{:});
     sd{m, n} = vertcat(SD{m, n}{:});
     mav{m, n} = vertcat(MAV{m, n}{:});
     iemg{m, n} = vertcat(IEMG{m, n}{:});

    disp(m)
    end
    disp(m)
end

%% Storing vertcat feature matrix's one row according to the subject
rms1 = rms(1,:);
rms2 = rms(2,:);
rms3 = rms(3,:);
rms4 = rms(4,:);
rms5 = rms(5,:);
sd1 = sd(1,:);
sd2 = sd(2,:);
sd3 = sd(3,:);
sd4 = sd(4,:);
sd5 = sd(5,:);
asm1 = asm(1,:);
asm2 = asm(2,:);
asm3 = asm(3,:);
asm4 = asm(4,:);
asm5 = asm(5,:);
ass1 = ass(1,:);
ass2 = ass(2,:);
ass3 = ass(3,:);
ass4 = ass(4,:);
ass5 = ass(5,:);
mav1 = mav(1,:);
mav2 = mav(2,:);
mav3 = mav(3,:);
mav4 = mav(4,:);
mav5 = mav(5,:);
mad1 = mad(1,:);
mad2 = mad(2,:);
mad3 = mad(3,:);
mad4 = mad(4,:);
mad5 = mad(5,:);
iemg1 = iemg(1,:);
iemg2 = iemg(2,:);
iemg3 = iemg(3,:);
iemg4 = iemg(4,:);
iemg5 = iemg(5,:);

%% Subject wise vertcat features stored in the form (subject1 - t1)
t1ass = vertcat(ass1{:});
t2ass = vertcat(ass2{:});
t3ass = vertcat(ass3{:});
t4ass = vertcat(ass4{:});
t5ass= vertcat(ass5{:});
t1mav = vertcat(mav1{:});
t2mav = vertcat(mav2{:});
t3mav = vertcat(mav3{:});
t4mav = vertcat(mav4{:});
t5mav = vertcat(mav5{:});
t1mad = vertcat(mad1{:});
t2mad = vertcat(mad2{:});
t3mad = vertcat(mad3{:});
t4mad = vertcat(mad4{:});
t5mad = vertcat(mad5{:});
t1iemg = vertcat(iemg1{:});
t2iemg = vertcat(iemg2{:});
t3iemg = vertcat(iemg3{:});
t4iemg = vertcat(iemg4{:});
t5iemg = vertcat(iemg5{:});
t1sd = vertcat(sd1{:});
t2sd = vertcat(sd2{:});
t3sd = vertcat(sd3{:});
t4sd = vertcat(sd4{:});
t5sd = vertcat(sd5{:});
t1rms = vertcat(rms1{:});
t2rms = vertcat(rms2{:});
t3rms = vertcat(rms3{:});
t4rms = vertcat(rms4{:});
t5rms = vertcat(rms5{:});
t1asm = vertcat(asm1{:});
t2asm = vertcat(asm2{:});
t3asm = vertcat(asm3{:});
t4asm = vertcat(asm4{:});
t5asm = vertcat(asm5{:});

%% Labelling of 2 types(binary, 14 class)
%v - subject one with 14 class classification
%ppv - subject one with power precision labelling
v = [ones(251,1);3*ones(246,1);4*ones(252,1);5*ones(246,1);6*ones(246,1);10*ones(252,1);20*ones(247,1);7*ones(252,1);12*ones(252,1);13*ones(252,1);15*ones(247,1);16*ones(246,1);18*ones(252,1);9*ones(252,1)];

w = [ones(253,1);3*ones(246,1);4*ones(252,1);5*ones(248,1);6*ones(249,1);10*ones(252,1);20*ones(246,1);7*ones(253,1);12*ones(252,1);13*ones(252,1);15*ones(247,1);16*ones(247,1);18*ones(252,1);9*ones(252,1)];

x = [ones(250,1);3*ones(246,1);4*ones(252,1);5*ones(247,1);6*ones(247,1);10*ones(252,1);20*ones(247,1);7*ones(252,1);12*ones(252,1);13*ones(252,1);15*ones(246,1);16*ones(247,1);18*ones(252,1);9*ones(251,1)];

y = [ones(252,1);3*ones(248,1);4*ones(253,1);5*ones(247,1);6*ones(247,1);10*ones(252,1);20*ones(246,1);7*ones(252,1);12*ones(252,1);13*ones(252,1);15*ones(246,1);16*ones(248,1);18*ones(252,1);9*ones(252,1)];

z = [ones(249,1);3*ones(248,1);4*ones(252,1);5*ones(247,1);6*ones(246,1);10*ones(252,1);20*ones(246,1);7*ones(252,1);12*ones(252,1);13*ones(252,1);15*ones(246,1);16*ones(247,1);18*ones(252,1);9*ones(252,1)];

ppv = [zeros(1740,1) ; ones(1753,1)];
ppw = [zeros(1746,1) ; ones(1755,1)];
ppx = [zeros(1741,1) ; ones(1752,1)];
ppy = [zeros(1745,1) ; ones(1754,1)];
ppz = [zeros(1740,1) ; ones(1753,1)];

% then we horizontally concatenated following labels with subject features accordingly
 
%% combining 7 features horizontally for a particular subject
%following is for subject 5
%ppl - power precision labelling
%l - 14 class classification
x1 = t5asm;
x2 = t5ass;
x3 = t5rms;
x4 = t5mad;
x5 = t5sd;
x6= t5mav;
x7= t5iemg;
combined5 = cell2mat([x1,x2,x3,x4,x5,x6,x7,]);
l_combined5 = [combined5,z];
ppl_combined5 =  [combined5,ppz];
        
