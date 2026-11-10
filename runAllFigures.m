function runAllFigures()
% Run every figure script and save each figure to results/ as .png and .fig.
% Usage: cd to this folder, then type  runAllFigures

scripts = {'secondFigure', 'thirdFigure', 'fourthFigure', 'fifthFigure', ...
           'sixthFigure', 'seventhFigure', 'eighthFigure', 'ninthFigure', ...
           'tenthFigure', 'eleventhFigure', 'twelfthFigure'};
names   = {'2nd figure', '3rd figure', '4th figure', '5th figure', ...
           '6th figure', '7th figure', '8th figure', '9th figure', ...
           '10th figure', '11th figure', '12th figure'};

outdir = fullfile(fileparts(mfilename('fullpath')), 'results');
if ~exist(outdir, 'dir'), mkdir(outdir); end

close all;
failed = {};
tAll = tic;
for k = 1:numel(scripts)
    fprintf('%-15s ... ', scripts{k});
    t = tic;
    try
        runOne(scripts{k});
        fig = gcf;
        try fig.Theme = 'light'; catch, end        % white background in saved files (R2025a+)
        png = fullfile(outdir, [names{k} '.png']);
        if exist('exportgraphics', 'file')
            exportgraphics(fig, png, 'Resolution', 200);
        else
            saveas(fig, png);
        end
        savefig(fig, fullfile(outdir, [names{k} '.fig']));
        close(fig);
        fprintf('saved (%.1f s)\n', toc(t));
    catch err
        fprintf('FAILED: %s\n', err.message);
        failed{end+1} = scripts{k}; %#ok<AGROW>
    end
end

fprintf('\nDone in %.1f s. Output folder: %s\n', toc(tAll), outdir);
if ~isempty(failed)
    fprintf('Failed: %s\n', strjoin(failed, ', '));
end
end

function runOne(name)
% Run a script inside its own function workspace, so the "clc; clear;"
% at the top of each script cannot wipe the loop variables above.
run(name);
end
