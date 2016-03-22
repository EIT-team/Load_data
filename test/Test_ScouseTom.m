
% datadir='C:\Users\James\Documents\Test Data\';
datadir='E:\Load_data_testing\';

%% load bdfs

somethingwentwrong=0;


try
    BVs=ScouseTom_LoadBV([datadir 'Baseline1.bdf']);
catch err
    fprintf(2,'error on initial test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'Baseline1_missingstart.bdf']);
catch err
    fprintf(2,'error on missing start test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'Baseline1_missingend.bdf']);
catch err
    fprintf(2,'error on missing end test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'Baseline1_missingboth.bdf']);
catch err
    fprintf(2,'error on missing both test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');
%% multifreq

try
    BVs=ScouseTom_LoadBV([datadir 'MultiFreq.bdf']);
catch err
    fprintf(2,'error on Multifreq test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'MultiFreq_missingstart.bdf']);
catch err
    fprintf(2,'error on Multifreq missing start test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'MultiFreq_missingend.bdf']);
catch err
    fprintf(2,'error on Multifreq missing end test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'MultiFreq_missingboth.bdf']);
catch err
    fprintf(2,'error on Multifreq missing both test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

%% Z check
try
    Z=ScouseTom_LoadZ([datadir 'firstinnir_zcheck_1.bdf']);
catch err
    fprintf(2,'error on Z load \n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end


%% AC check

try
    BVs=ScouseTom_LoadBV([datadir 'rptest2.eeg']);
catch err
    fprintf(2,'error on Multifreq missing both test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_LoadBV([datadir 'rptest3.eeg']);
catch err
    fprintf(2,'error on Multifreq missing both test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');


%% Big file

try
    BVs=ScouseTom_LoadBV([datadir 'RealThing_NIRFACE_2.bdf']);
catch err
    fprintf(2,'error on Multifreq big file test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');



%% did it work?

if somethingwentwrong
    fprintf(2,'FUCK! SOMETHING WENT WRONG!\n');
else
    fprintf('YAY! IT WORKED :)\n');
end

