% test code is functioning on real examples. Requires files from https://doi.org/10.5281/zenodo.4783547

datadir='I:/Load_data_testing/';
%datadir='data/';
%% load bdfs

somethingwentwrong=0;


try
    BVs=ScouseTom_Load([datadir 'Baseline1.bdf']);
catch err
    fprintf(2,'error on initial test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_Load([datadir 'Baseline1_missingstart.bdf']);
catch err
    fprintf(2,'error on missing start test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_Load([datadir 'Baseline1_missingend.bdf']);
catch err
    fprintf(2,'error on missing end test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_Load([datadir 'Baseline1_missingboth.bdf']);
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
    BVs=ScouseTom_Load([datadir 'MultiFreq.bdf']);
catch err
    fprintf(2,'error on Multifreq test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_Load([datadir 'MultiFreq_missingstart.bdf']);
catch err
    fprintf(2,'error on Multifreq missing start test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_Load([datadir 'MultiFreq_missingend.bdf']);
catch err
    fprintf(2,'error on Multifreq missing end test\n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_Load([datadir 'MultiFreq_missingboth.bdf']);
catch err
    fprintf(2,'error on Multifreq missing both test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

%% Z check
try
    Z=ScouseTom_Load([datadir 'firstinnir_zcheck_1.bdf']);
catch err
    fprintf(2,'error on Z load \n');
    somethingwentwrong=1;
    fprintf(2, '%s\n', getReport(err, 'extended'));
end


%% AC check

try
    BVs=ScouseTom_Load([datadir 'rptest2.eeg']);
catch err
    fprintf(2,'error on Actichamp test 1\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');

try
    BVs=ScouseTom_Load([datadir 'rptest3.eeg']);
catch err
    fprintf(2,'error on Actichamp test 2\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');


%% Big file

try
    BVs=ScouseTom_Load([datadir 'RealThing_NIRFACE_2.bdf']);
catch err
    fprintf(2,'error on Multifreq big file test\n');
    somethingwentwrong=1;
end

disp('=================================================');
disp('=================================================');
disp('=================================================');



%% did it work?

if somethingwentwrong
    fprintf(2,'Oh No! SOMETHING WENT WRONG!\n');
else
    fprintf('YAY! IT WORKED :)\n');
end
