classdef CMEApp < handle
%CMEAPP  Interactive GUI for the BME 404 CME segmentation project.
%
% USAGE
%   >> cd('<project root>'); addpath('src','gui');  app = CMEApp;
%   or simply   >> cmeGUI
%
% WHAT IT INTEGRATES (the course requires the GUI to cover every feature)
%   * load any B-scan of any Duke subject (dropdown + slider)
%   * switch between the published BASELINE pipeline and OUR MODIFIED one
%   * run the full segmentation and show the contour overlaid on the B-scan
%   * step through every intermediate stage of the pipeline for explanation
%     during the viva (denoised, enhanced, ROI, each of the four directional
%     contour maps, the fused map, the candidate basins, the final mask)
%   * show the clinical read-out: fluid area, cyst count, mean diameter,
%     central subfield thickness, and the severity grade
%   * show Dice / precision / recall against both expert graders whenever the
%     selected B-scan is one of the 78 with expert annotation
%   * export the current view as a PNG
%
% ARCHITECTURE NOTE
%   This class contains NO algorithm code.  Every callback calls the same
%   functions in src/ that the batch scripts use (segmentCME, evalSegmentation,
%   clinicalQuantify, ...).  Nothing is duplicated, so the GUI can never drift
%   out of step with the results in the report.
%
% Written as a programmatic uifigure rather than a binary .mlapp so that the
% whole interface is readable, reviewable, diff-able source code -- which also
% makes it far easier to explain in the viva.

    properties (Access = public)
        Fig
        SubjectDD
        BscanSlider
        BscanLabel
        ModeSwitch
        StageDD
        GradedDD
        RunBtn
        ExportBtn
        AutoChk
        AxRaw
        AxMain
        ClinicalArea
        MetricArea
        StatusLbl
        InfoArea
    end

    properties (Access = private)
        Scan        = []     % current loadBScan struct
        Result      = []     % current segmentCME output
        EvalIndex   = []     % buildEvalSet table
        Busy        = false
    end

    methods
        function app = CMEApp()
            app.ensurePath();
            app.EvalIndex = buildEvalSet();
            app.buildUI();
            app.refreshGradedList();
            app.loadCurrentScan();
            app.runPipeline();
        end

        function delete(app)
            if isvalid(app.Fig), delete(app.Fig); end
        end
    end

    % =====================================================================
    % These are public rather than private on purpose: it lets the pipeline be
    % driven from the command line during the demo (app.runPipeline(), etc.)
    % and lets the GUI be exercised by a test script.
    methods (Access = public)

        function ensurePath(~)
            here = fileparts(mfilename('fullpath'));       % ...\gui
            root = fileparts(here);
            if exist(fullfile(root,'src'),'dir')
                addpath(fullfile(root,'src'));
            end
        end

        % -----------------------------------------------------------------
        function buildUI(app)
            app.Fig = uifigure('Name','BME 404 - CME Segmentation and Severity Grading (A2 Group 03)', ...
                'Position',[60 60 1500 900], 'Color',[0.96 0.96 0.97]);

            G = uigridlayout(app.Fig, [1 2]);
            G.ColumnWidth = {320, '1x'};
            G.RowHeight   = {'1x'};

            % ---------- left control column ------------------------------
            left = uigridlayout(G,[7 1]);
            % Panel 1 holds three controls (subject, B-scan slider, jump-to
            % list) and a slider needs vertical room for its tick labels, so
            % it gets noticeably more height than the others.
            left.RowHeight = {150, 115, 78, 155, 120, '1x', 22};
            left.Padding = [8 8 8 8];

            % -- data selection
            pData = uipanel(left,'Title','1. Data','FontWeight','bold');
            gd = uigridlayout(pData,[3 2]);
            gd.ColumnWidth = {70,'1x'};
            gd.RowHeight   = {22, 46, 22};   % the slider row needs the space
            gd.RowSpacing  = 6;
            uilabel(gd,'Text','Subject');
            app.SubjectDD = uidropdown(gd, ...
                'Items', arrayfun(@(k) sprintf('Subject_%02d',k), 1:10, 'UniformOutput', false), ...
                'ItemsData', 1:10, 'Value', 5, ...
                'ValueChangedFcn', @(s,e) app.onSubjectChanged());
            app.BscanLabel = uilabel(gd,'Text','B-scan 24');
            app.BscanSlider = uislider(gd,'Limits',[1 61],'Value',24, ...
                'MajorTicks',[1 15 30 45 61], ...
                'ValueChangedFcn', @(s,e) app.onScanChanged());
            uilabel(gd,'Text','Jump to');
            % Only 110 of the 610 B-scans were graded.  Stepping the slider at
            % random almost always lands on an ungraded scan, which makes the
            % tool look wrong.  This dropdown goes straight to the ones that
            % have ground truth.
            app.GradedDD = uidropdown(gd, 'Items', {'(any B-scan)'}, ...
                'ValueChangedFcn', @(s,e) app.onGradedPick());

            % -- pipeline mode
            pMode = uipanel(left,'Title','2. Pipeline','FontWeight','bold');
            gm = uigridlayout(pMode,[3 2]); gm.ColumnWidth = {110,'1x'};
            uilabel(gm,'Text','Mode');
            app.ModeSwitch = uidropdown(gm, ...
                'Items',{'modified (ours)','baseline (Liu 2021)'}, ...
                'ItemsData',{'modified','baseline'}, 'Value','modified', ...
                'ValueChangedFcn', @(s,e) app.runPipeline());
            uilabel(gm,'Text','Auto-run');
            app.AutoChk = uicheckbox(gm,'Text','on scan change','Value',true);
            app.RunBtn = uibutton(gm,'Text','Run segmentation', ...
                'BackgroundColor',[0.20 0.45 0.75],'FontColor','w','FontWeight','bold', ...
                'ButtonPushedFcn', @(s,e) app.runPipeline());
            uibutton(gm,'Text','Export view', ...
                'ButtonPushedFcn', @(s,e) app.exportView());

            % -- stage viewer
            pStage = uipanel(left,'Title','3. Intermediate stage','FontWeight','bold');
            gs = uigridlayout(pStage,[1 1]);
            app.StageDD = uidropdown(gs, 'Items', app.stageNames(), 'Value','final overlay', ...
                'ValueChangedFcn', @(s,e) app.refreshMain());

            % -- clinical output
            pClin = uipanel(left,'Title','4. Clinical output (our addition)','FontWeight','bold');
            gc = uigridlayout(pClin,[1 1]);
            app.ClinicalArea = uitextarea(gc,'Editable','off','FontName','Consolas','FontSize',11);

            % -- metrics
            pMet = uipanel(left,'Title','5. Accuracy vs expert graders','FontWeight','bold');
            gme = uigridlayout(pMet,[1 1]);
            app.MetricArea = uitextarea(gme,'Editable','off','FontName','Consolas','FontSize',11);

            % -- explanation of the current stage
            pInfo = uipanel(left,'Title','6. What you are looking at','FontWeight','bold');
            gi = uigridlayout(pInfo,[1 1]);
            app.InfoArea = uitextarea(gi,'Editable','off','FontSize',10.5,'WordWrap','on');

            app.StatusLbl = uilabel(left,'Text','ready','FontAngle','italic');

            % ---------- right image column -------------------------------
            right = uigridlayout(G,[2 1]);
            right.RowHeight = {'1x','1.4x'};

            pRaw = uipanel(right,'Title','Raw B-scan','FontWeight','bold');
            gr = uigridlayout(pRaw,[1 1]);
            app.AxRaw = uiaxes(gr); app.AxRaw.Toolbar.Visible = 'off';

            pMain = uipanel(right,'Title','Result','FontWeight','bold');
            gmn = uigridlayout(pMain,[1 1]);
            app.AxMain = uiaxes(gmn);
        end

        % -----------------------------------------------------------------
        function n = stageNames(~)
            n = {'final overlay', ...
                 'raw', 'denoised', 'contrast enhanced', ...
                 'ROI (retina + inner band)', 'row-mean-gray curve', ...
                 'wave potential energy', 'wave map (sea/wave/coast)', ...
                 'direction 1  (theta = 0)', 'direction 2  (theta = pi/2)', ...
                 'direction 3  (theta = pi)', 'direction 4  (theta = 3pi/2)', ...
                 'fused contours', 'candidate basins', 'final mask only', ...
                 'expert A vs expert B'};
        end

        % -----------------------------------------------------------------
        function refreshGradedList(app)
            % Build the list of graded B-scans for the selected subject.
            T = app.EvalIndex;
            rows = T(T.subject == app.SubjectDD.Value & logical(T.hasLayers), :);
            items = {'(any B-scan)'};
            data  = {[]};
            for i = 1:height(rows)
                if rows.nFluid1(i) > 0
                    items{end+1} = sprintf('k%02d  fluid %d px', rows.bscan(i), rows.nFluid1(i)); %#ok<AGROW>
                else
                    items{end+1} = sprintf('k%02d  graded, no fluid', rows.bscan(i)); %#ok<AGROW>
                end
                data{end+1} = rows.bscan(i); %#ok<AGROW>
            end
            app.GradedDD.Items     = items;
            app.GradedDD.ItemsData = data;
            app.GradedDD.Value     = data{1};
        end

        function onSubjectChanged(app)
            app.refreshGradedList();
            app.onScanChanged();
        end

        function onGradedPick(app)
            v = app.GradedDD.Value;
            if isempty(v), return; end
            app.BscanSlider.Value = v;
            app.onScanChanged();
        end

        function onScanChanged(app)
            app.BscanSlider.Value = round(app.BscanSlider.Value);
            app.BscanLabel.Text = sprintf('B-scan %d', app.BscanSlider.Value);
            app.loadCurrentScan();
            if app.AutoChk.Value
                app.runPipeline();
            else
                app.refreshMain();
            end
        end

        function loadCurrentScan(app)
            s = app.SubjectDD.Value;
            k = round(app.BscanSlider.Value);
            app.setStatus(sprintf('loading Subject_%02d B-scan %d ...', s, k));
            try
                app.Scan = loadBScan(s, k);
            catch ME
                uialert(app.Fig, ME.message, 'Could not load data');
                app.setStatus('load failed');
                return
            end
            imshow(app.Scan.img, [0 255], 'Parent', app.AxRaw);
            if app.Scan.isGraded
                gtxt = sprintf('GRADED - expert A marked %d fluid px', nnz(app.Scan.fluid1));
            else
                gtxt = 'NOT GRADED by the experts (no ground truth exists here)';
            end
            title(app.AxRaw, sprintf('Subject_%02d, B-scan %d   [%s]', s, k, gtxt));
            app.setStatus('loaded');
        end

        % -----------------------------------------------------------------
        function runPipeline(app)
            if app.Busy || isempty(app.Scan), return; end
            app.Busy = true;
            app.RunBtn.Enable = 'off';
            app.setStatus('running the pipeline ...'); drawnow;
            try
                cfg = cmeConfig(app.ModeSwitch.Value);
                app.Result = segmentCME(app.Scan.img, cfg, app.Scan.stack);
                app.updateClinical();
                app.updateMetrics();
                app.refreshMain();
                app.setStatus(sprintf('done in %.2f s  (%s mode)', ...
                    app.Result.timing.total, app.ModeSwitch.Value));
            catch ME
                uialert(app.Fig, sprintf('%s\n\n(%s line %d)', ME.message, ...
                    ME.stack(1).name, ME.stack(1).line), 'Pipeline error');
                app.setStatus('error');
            end
            app.RunBtn.Enable = 'on';
            app.Busy = false;
        end

        % -----------------------------------------------------------------
        function updateClinical(app)
            C = app.Result.clinical;
            app.ClinicalArea.Value = { ...
                sprintf('fluid area   : %7.0f px', C.areaPx), ...
                sprintf('             : %7.4f mm^2', C.areaMm2), ...
                sprintf('cyst count   : %7d', C.cystCount), ...
                sprintf('mean diameter: %7.0f um', C.meanDiamUm), ...
                sprintf('max  diameter: %7.0f um', C.maxDiamUm), ...
                sprintf('CST          : %7.0f um', C.cstUm), ...
                sprintf('retina thick.: %7.0f um', C.retinalThicknessUm), ...
                '', ...
                sprintf('GRADE (area) : %s', upper(C.gradeDetail.areaGrade)), ...
                sprintf('GRADE (CST)  : %s', upper(C.gradeDetail.cstGrade)), ...
                sprintf('reported     : %s', upper(C.grade))};
        end

        function updateMetrics(app)
            S = app.Scan;

            % The Duke experts graded only 110 of the 610 B-scans.  On the
            % other 500 the fluid mask is entirely NaN -- those scans were
            % never looked at, so "we segmented something and the expert did
            % not" is meaningless there.  Saying so explicitly is important:
            % otherwise scrubbing through the volume looks full of false
            % positives when most of it simply has no ground truth.
            if ~S.isGraded
                app.MetricArea.Value = { ...
                    'THIS B-SCAN WAS NEVER GRADED.', ...
                    '', ...
                    'Only 110 of the 610 Duke B-scans carry', ...
                    'expert annotation (11 per subject); on the', ...
                    'rest manualFluid1 is entirely NaN.', ...
                    '', ...
                    'No score can be computed here, and our', ...
                    'output is NOT a false positive -- there is', ...
                    'simply nothing to compare it against.', ...
                    '', ...
                    sprintf('our segmented area: %d px', nnz(app.Result.mask))};
                return
            end

            if ~any(S.fluid1(:)) && ~any(S.fluid2(:))
                app.MetricArea.Value = { ...
                    'GRADED, and both experts found NO fluid.', ...
                    '', ...
                    sprintf('our false-positive area: %d px', nnz(app.Result.mask)), ...
                    '', ...
                    'This is a genuine negative: a scan the', ...
                    'experts examined and declared fluid-free.', ...
                    'Only 32 B-scans are in this category.'};
                return
            end
            valid = app.Result.roi.mask & S.gradedMask;
            pred  = app.Result.mask & valid;
            EA = evalSegmentation(pred, S.fluid1, valid);
            EB = evalSegmentation(pred, S.fluid2, valid);
            EX = evalSegmentation(S.fluid2, S.fluid1, valid);
            app.MetricArea.Value = { ...
                '              Dice   Prec  Recall', ...
                sprintf('vs expert A : %.3f  %.3f  %.3f', EA.dice, EA.precision, EA.recall), ...
                sprintf('vs expert B : %.3f  %.3f  %.3f', EB.dice, EB.precision, EB.recall), ...
                '', ...
                sprintf('A vs B (inter-observer): %.3f', EX.dice), ...
                '(that is the practical ceiling)'};
        end

        % -----------------------------------------------------------------
        function refreshMain(app)
            if isempty(app.Result), return; end
            R = app.Result; S = app.Scan;
            ax = app.AxMain; cla(ax, 'reset');

            % crop to the retina so detail is visible
            r0 = max(1, round(min(R.roi.ilm)) - 30);
            r1 = min(size(S.img,1), round(max(R.roi.obm)) + 30);
            cr = @(X) X(r0:r1, :, :);
            valid = R.roi.mask & S.gradedMask;

            stage = app.StageDD.Value;
            info  = '';

            switch stage
                case 'final overlay'
                    imshow(vizOverlay(cr(S.img), {cr(S.fluid1), cr(S.fluid2), cr(R.mask & valid)}, ...
                        {[0 1 0],[0.2 0.6 1],[1 0 0]}, struct('fillAlpha',0.20,'lineWidth',2)),'Parent',ax);
                    title(ax,'green = expert A, blue = expert B, red = our segmentation');
                    info = ['The final CME contour on the original B-scan, with both expert ' ...
                            'tracings for comparison. Filled red is what the pipeline reports ' ...
                            'as fluid.'];

                case 'raw'
                    imshow(cr(R.stages.raw),[0 255],'Parent',ax); title(ax,'raw B-scan');
                    info = 'The unprocessed OCT B-scan as stored in the Duke dataset.';

                case 'denoised'
                    imshow(cr(R.stages.denoised),[0 255],'Parent',ax);
                    title(ax, sprintf('denoised (%s)', R.cfg.denoise.method));
                    info = ['Speckle suppression. The modified pipeline uses a 5x5 MEDIAN ' ...
                            'filter, which is a rank filter: it discards outlier speckle ' ...
                            'spikes instead of averaging them in, so step edges survive. ' ...
                            'The baseline uses a Gaussian blur, which is linear and softens ' ...
                            'the very edges the wave operator needs.'];

                case 'contrast enhanced'
                    imshow(cr(R.stages.enhanced),[0 255],'Parent',ax);
                    title(ax, sprintf('contrast enhanced (%s)', R.cfg.contrast.method));
                    info = ['CLAHE equalises the histogram of each tile separately and ' ...
                            'interpolates between tiles, so it adapts to the uneven ' ...
                            'illumination of an OCT scan. The clip limit stops it amplifying ' ...
                            'speckle. Tile count matters: a tile must be larger than a cyst.'];

                case 'ROI (retina + inner band)'
                    im = vizOverlay(cr(R.work), {cr(R.roi.mask), cr(R.roi.innerBand)}, ...
                        {[0 0.8 1],[1 0.9 0]}, struct('fillAlpha',0.18));
                    imshow(im,'Parent',ax); hold(ax,'on');
                    plot(ax, R.roi.ilm - r0 + 1, 'g','LineWidth',1.5);
                    plot(ax, R.roi.obm - r0 + 1, 'y','LineWidth',1.5);
                    hold(ax,'off');
                    title(ax, sprintf('ROI: retina %.0f px thick, inner band %.2f-%.2f', ...
                        R.roi.thicknessPx, R.cfg.roi.depthBand));
                    info = ['Cyan is the retina between ILM and OBM; yellow is the inner ' ...
                            'band where CME actually occurs. Restricting the search here is ' ...
                            'the "position feature" screening of the paper, and it removes ' ...
                            'the dark choroidal vessels that would otherwise look like cysts.'];

                case 'row-mean-gray curve'
                    plot(ax, R.roi.rowMeanRaw, 1:numel(R.roi.rowMeanRaw), 'Color',[.7 .7 .7]);
                    hold(ax,'on');
                    plot(ax, R.roi.rowMean, 1:numel(R.roi.rowMean),'k','LineWidth',1.5);
                    xline(ax, R.roi.imgMean,'b','LineWidth',2);
                    yline(ax, R.roi.A,'r--','LineWidth',1.5);
                    yline(ax, R.roi.B,'r--','LineWidth',1.5);
                    set(ax,'YDir','reverse'); grid(ax,'on'); hold(ax,'off');
                    xlabel(ax,'mean gray of row'); ylabel(ax,'row');
                    title(ax, sprintf('Eq.(1): A = %d, B = %d', R.roi.A, R.roi.B));
                    info = ['Equation (1) of the paper. The mean gray of every row is ' ...
                            'plotted against depth; where it crosses the whole-image mean ' ...
                            '(blue) gives points A and B, which bracket the retina.'];

                case 'wave potential energy'
                    imagesc(ax, cr(R.dirs.exemplar.wave)); axis(ax,'image'); axis(ax,'off');
                    colorbar(ax); title(ax,'wave = \phi_v + \phi_g   (\theta = 0)');
                    info = ['The wave potential energy map. phi_g is the Gaussian-weighted ' ...
                            'gray value (gravitational term) and phi_v = v*vq*sigma is the ' ...
                            'kinetic term, which peaks where the gray value ramps steadily ' ...
                            'in the sweep direction.'];

                case 'wave map (sea/wave/coast)'
                    W = R.dirs.exemplar;
                    wm = zeros([numel(r0:r1) size(S.img,2) 3]);
                    wm(:,:,3) = double(cr(~W.IT & ~W.Pen));
                    wm(:,:,2) = double(cr(W.IT));
                    wm(:,:,1) = double(cr(W.Pen & ~W.IT));
                    imshow(wm,'Parent',ax);
                    title(ax,'blue = seawater, green = wave area, red = coastline');
                    info = ['The three-region picture from the 2020 paper. Sweeping forward ' ...
                            'the operator crosses seawater (low energy), the wave (high ' ...
                            'energy just before a boundary) and the coastline (inside the ' ...
                            'object). The boundary is the junction between wave and coast.'];

                case {'direction 1  (theta = 0)','direction 2  (theta = pi/2)', ...
                      'direction 3  (theta = pi)','direction 4  (theta = 3pi/2)'}
                    t = str2double(stage(11));
                    if t > numel(R.dirs.dirBnd), t = 1; end
                    imshow(vizOverlay(cr(R.work), {cr(R.dirs.dirBnd{t})}, {[1 0.25 0.25]}, ...
                        struct('lineWidth',1)),'Parent',ax);
                    title(ax, sprintf('contour found sweeping at %d degrees  (%.2f%% of ROI)', ...
                        R.dirs.anglesDeg(t), 100*mean(R.dirs.dirBnd{t}(R.roi.mask))));
                    info = ['One operating direction on its own. The operator only fires ' ...
                            'where the image goes dark-to-bright ALONG the sweep, so a ' ...
                            'downward sweep finds cyst floors, an upward sweep finds roofs, ' ...
                            'and the sideways sweeps find the walls. No single direction ' ...
                            'closes the contour, which is the whole argument of the paper.'];

                case 'fused contours'
                    imshow(vizOverlay(cr(R.work), {cr(R.dirs.fused)}, {[1 0.9 0.1]}, ...
                        struct('lineWidth',1)),'Parent',ax);
                    title(ax,'union of all four directions');
                    info = ['The four directional results superimposed: a pixel marked as a ' ...
                            'boundary in ANY direction becomes part of the omnidirectional ' ...
                            'contour.'];

                case 'candidate basins'
                    imshow(vizOverlay(cr(R.work), {cr(R.integrate.candidates)}, {[0 1 1]}, ...
                        struct('fillAlpha',0.35,'lineWidth',1)),'Parent',ax);
                    title(ax,'candidate regions before screening');
                    info = ['Candidate cyst regions before the gray, size, width and ' ...
                            'contrast screening are applied. Compare with the final mask to ' ...
                            'see what the screening removed.'];

                case 'final mask only'
                    imshow(cr(R.mask & valid),'Parent',ax);
                    title(ax, sprintf('final binary mask (%d px)', nnz(R.mask & valid)));
                    info = 'The binary segmentation the pipeline outputs.';

                case 'expert A vs expert B'
                    imshow(vizOverlay(cr(S.img), {cr(S.fluid1), cr(S.fluid2)}, ...
                        {[0 1 0],[0.2 0.6 1]}, struct('fillAlpha',0.20,'lineWidth',2)),'Parent',ax);
                    EX = evalSegmentation(S.fluid2, S.fluid1, valid);
                    title(ax, sprintf('the two graders agree at Dice %.3f', EX.dice));
                    info = ['The two expert graders compared with each other. Their ' ...
                            'disagreement is the realistic ceiling for any automatic method ' ...
                            'scored against one of them.'];
            end
            app.InfoArea.Value = info;
        end

        % -----------------------------------------------------------------
        function exportView(app)
            p = dukePaths();
            fn = fullfile(p.figures, sprintf('gui_export_S%02d_k%02d_%s.png', ...
                app.SubjectDD.Value, round(app.BscanSlider.Value), ...
                matlab.lang.makeValidName(app.StageDD.Value)));
            exportgraphics(app.AxMain, fn, 'Resolution', 200);
            app.setStatus(sprintf('exported %s', fn));
        end

        function setStatus(app, msg)
            app.StatusLbl.Text = msg;
            drawnow limitrate;
        end
    end
end
