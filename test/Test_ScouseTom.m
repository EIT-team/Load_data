
datadir='C:\Users\James\Documents\Test Data\';
% datadir='E:\testperchn\';

%% load bdfs

somethingwentwrong=0;


try
    BVs=ScouseTom_LoadBV([datadir 'Baseline1.bdf']);
catch
    fprintf(2,'error on initial test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'Baseline1_missingstart.bdf']);
catch
    fprintf(2,'error on missing start test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'Baseline1_missingend.bdf']);
catch
    fprintf(2,'error on missing end test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'Baseline1_missingboth.bdf']);
catch
    fprintf(2,'error on missing both test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');
%% multifreq

try
    BVs=ScouseTom_LoadBV([datadir 'MultiFreq.bdf']);
catch
    fprintf(2,'error on Multifreq test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'MultiFreq_missingstart.bdf']);
catch
    fprintf(2,'error on Multifreq missing start test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'MultiFreq_missingend.bdf']);
catch
    fprintf(2,'error on Multifreq missing end test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'MultiFreq_missingboth.bdf']);
catch
    fprintf(2,'error on Multifreq missing both test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

%% Z check
try
    Z=ScouseTom_LoadZ([datadir 'firstinnir_zcheck_1.bdf']);
catch
    fprintf(2,'error on Z load \n');
    somethingwentwrong=1;
end

%% did it work?

if somethingwentwrong
    fprintf(2,'FUCK! SOMETHING WENT WRONG!\n');
else
    fprintf('YAY! IT WORKED :)');
end

