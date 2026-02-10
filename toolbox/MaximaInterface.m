classdef MaximaInterface < handle
    % MaximaInterface
    % Persistent interface to a Maxima console instance
    % Cross-platform: Windows & Linux
    %
    % Example:
    %   maxima = MaximaInterface();
    %   result = maxima.sendAndParse('factor(x^2 - 1);');
    %   maxima.sendNoWait('a: 2$');
    %   res2 = maxima.sendAndParse('a + b;');

    properties
        maxWait = 10.0; % seconds
        startedAt = 0
        transcript_file = ''
    end

    properties (SetAccess = private)
        startupLines
        transcript_fid
    end

    properties (Access = private)
        proc        % Java Process object
        stdin       % OutputStreamWriter for Maxima stdin
        stdout      % BufferedReader for Maxima stdout
        nextId % Last used %i/%o number
    end

    properties (Access = private, Constant)
        PROMPT_INPUT = '(%i';
        PROMPT_OUTPUT = '(%o';
    end

    methods
        function obj = MaximaInterface(maxWait, transcript_file)
            arguments
                maxWait (1,1) double {mustBePositive} = 2.0
                transcript_file = ''
            end
            
            obj.nextId = 1;
            obj.maxWait = maxWait;
            obj.transcript_file = transcript_file;
            if ~isempty(transcript_file)
                obj.transcript_fid = fopen(transcript_file, 'w');
            end

            obj.startMaxima();
            obj.startedAt = datetime('now');

            obj.sendNoWait('load(lrats)');
            % obj.sendNoWait('load(gentran)');
            % obj.sendNoWait('gentranlang:''c$');
            % obj.sendNoWait('genfloat:true$');
            % obj.sendNoWait('clinelen:1000$');
            obj.sendNoWait('matchdeclare(exp_to_e_match, true)$');
            obj.sendNoWait('defrule(e_to_exp_rule, %e^x, exp(x))$');
            obj.sendNoWait('optimizeWithIntermediates(l, temp):= block([o, temps, pre_replace_list, replace_list, temp_exprs, temp_vars, i], o: optimize(l), if op(o)#''block then return([0, [], [], o]), temps: inpart(o, 1), pre_replace_list: makelist(temps[i]=concat(temp, i), i, 1, length(temps)), replace_list: [], temp_vars: [], temp_exprs: [], for i: 1 thru length(temps) do ( if not(listp(inpart(o, i+1, 2)) or matrixp(inpart(o, i+1, 2))) then ( temp_vars: endcons(concat(temp, i), temp_vars), temp_exprs: endcons(subst(pre_replace_list, inpart(o, i+1, 2)), temp_exprs), replace_list: endcons(pre_replace_list[i], replace_list) ) else ( replace_list: endcons(inpart(o, i+1, 1)=inpart(o, i+1, 2), replace_list) ) ), opts: subst(replace_list, subst(replace_list, inpart(o, length(temps)+2))), [length(temp_vars), temp_vars, temp_exprs, opts] )$');
        end

        function idStr = sendNoWait(obj, cmd)
            % Send a command to Maxima without output (terminated by $)
            % No further ; or $ must be in the command!
            % Returns the new output ID (%oN)
            if isempty(obj.proc) || ~obj.proc.isAlive()
                delete(obj);
                error('Maxima process has died. Interface and all expressions are now invalid.');
            end

            if ~endsWith(strtrim(cmd), '$')
                cmd = [cmd, '$'];
            end
            % % Safe mode
            % [~, extraLines, id] = sendAndParse(obj, cmd);
            % msg = strjoin(extraLines, newline);
            % if contains(msg, 'error') || contains(msg, 'incorrect')
            %     error('Maxima error: "%s"\nfrom command "%s".', msg, cmd)
            % end
            % idStr = sprintf('%%o%d', id);

            obj.stdin.write([cmd, newline]);
            obj.stdin.flush();

            % Maxima updates the internal %o number even without output
            idStr = sprintf('%%o%d', obj.nextId);
            obj.nextId = obj.nextId + 1;
        end

        function [result, extraLines, id] = sendAndParse(obj, cmd, expect_extra)
            arguments
                obj
                cmd
                expect_extra = false
            end
            % Send a command to Maxima and parse the response
            % Returns result string and output ID (%oN)
            if isempty(obj.proc) || ~obj.proc.isAlive()
                delete(obj);
                error('Maxima process has died. Interface and all expressions are now invalid.');
            end
    
            if ~isempty(cmd)
                id = obj.nextId;
                expect_outprompt = true;
                if endsWith(strtrim(cmd), '$')
                    expect_outprompt = false;
                elseif ~endsWith(strtrim(cmd), ';')
                    cmd = [cmd, ';'];
                end
                obj.stdin.write([cmd, newline]);
                obj.stdin.flush();
                
                obj.nextId = obj.nextId + 1;
            else
                id = [];
            end

            [result, extraLines] = obj.readToNextInputPrompt(expect_extra);

            if ~isempty(cmd) && expect_outprompt && isempty(result)
                warning('Was expecting output from Maxima but did not get any.')
            end
        end

        function delete(obj)
            % Cleanly stop Maxima when deleting the object
            if ~isempty(obj.proc) && obj.proc.isAlive()
                try
                    obj.stdin.write(['quit();', newline]);
                catch
                    % Ignore if process is already dead
                end
            end
            if ~isempty(obj.transcript_file)
                fclose(obj.transcript_fid);
            end
        end

        function checkAlive(obj)
            if isempty(obj.proc) || ~obj.proc.isAlive()
                error('MaximaExpression: The Maxima process has died. Expression is invalid.');
            end
        end            
        
    end

    methods (Static)
        function obj = getInstance(maxWait, transcript_file)
            arguments
                maxWait (1,1) double {mustBePositive} = 10.0
                transcript_file = ''
            end

            % Singleton accessor
            persistent instance
            if isempty(instance) || ~isvalid(instance)
                instance = MaximaInterface(maxWait, transcript_file);
            end
            obj = instance;
        end

        function stopInstance()
            % Stop and clear the singleton instance
            persistent instance
            if ~isempty(instance) && isvalid(instance)
                delete(instance);
            end
            instance = [];
        end
    end

    methods (Access = private)
        function [outputLines, extraLines] = readToNextInputPrompt(obj, expect_extra)
            arguments
                obj
                expect_extra = false
            end
            % Read available output until the all expected output is seen.
            % Collect lines without prompts (startup), output prompt lines,
            % discard empty lines and input prompt lines. Update nextId
            % based on the highest %o or %i seen.

            outputLines = {};
            extraLines = {};

            tStart = tic;
            while ~obj.stdout.ready() && toc(tStart) < obj.maxWait
                pause(0.01);
            end

            lastSeenIn = 0;
            continued_output = false;
            list_output = false;        % Lists don't have line continuation characters
            while toc(tStart) < obj.maxWait
                % Let's hope we always read entire lines
                raw = obj.readAll();
                if ~isempty(raw)
                    tStart = tic; % reset timeout on new data
                end
                lines = strsplit(raw, newline);
                for i = 1:length(lines)
                    line = strtrim(lines{i});
                    if isempty(line)
                        continue;
                    end
                    
                    % detect and process input prompt
                    [tok, id_end] = regexp(line, '^\(%i(\d+)\)', 'once', 'tokens', 'end');
                    if ~isempty(tok)
                        lastSeenIn = str2double(tok{1});
                        % We may have to catch up to the current input id
                        if lastSeenIn > obj.nextId
                            warning(['Maxima input prompt IDs are not sequential. Read %%i%d but current should be %%i%d.' newline ...
                                'Saw this output: "%s"'], lastSeenIn, obj.nextId, strjoin(extraLines, newline));
                            obj.nextId = lastSeenIn;
                        end
                        line = strtrim(line(id_end+1:end));
                        % continue;
                    end

                    % detect and process ouput prompt
                    [tok, id_end] = regexp(line, '^\(%o(\d+)\)', 'once', 'tokens', 'end');
                    if ~isempty(tok)
                        continued_output = endsWith(line, '\');
                        if continued_output
                            outputLine = strtrim(line(id_end+1:end-1));
                        else
                            outputLine = strtrim(line(id_end+1:end));
                        end
                        if startsWith(outputLine, '[') && ~endsWith(outputLine, ']')
                            list_output = true;
                        end
                        if isempty(outputLine)
                            continued_output = true; % some times the output is only in the next line
                        end

                        outputLines{end+1} = outputLine; %#ok<AGROW>
                        % Output ids may be more than one at a time. Next input id must be one more than last output id.
                        obj.nextId = str2double(tok{1}) + 1;
                    else
                        if continued_output
                            continued_output = endsWith(line, '\');
                            if continued_output % ignore '\' if it not at the end of a continuation line
                                outputLines{end} = [outputLines{end}, line(1:end-1)];
                            else
                                outputLines{end} = [outputLines{end}, line];
                            end
                        elseif list_output
                            outputLines{end} = [outputLines{end}, line];
                            list_output = ~endsWith(line, ']');
                        else
                            extraLines{end+1} = line; %#ok<AGROW>
                            expect_extra = false;
                        end
                    end
                end
                if lastSeenIn == obj.nextId && ~expect_extra
                    break;
                end
                pause(0.01);
            end
            if toc(tStart) >= obj.maxWait
                warning('Timeout while reading from Maxima process. Expression IDs may be out of sync.');
            end
            if continued_output || list_output
                warning('Continued output from Maxima (backslash or list) was not properly finished. This may cause parsing issues. Output was: "%s"', strjoin(outputLines, newline));
            end
            obj.nextId = lastSeenIn;
        end

        function startMaxima(obj)
            % Start a new Maxima instance
            import java.io.*
            import java.lang.*

            try
                % Start process
                if ispc
                    cmd = {'cmd.exe', '/c', 'maxima'};
                else
                    cmd = {'maxima'};
                end
                pb = ProcessBuilder(cmd);
                pb.redirectErrorStream(true);
                obj.proc = pb.start();

                % Initialize streams
                obj.stdin = BufferedWriter(OutputStreamWriter(obj.proc.getOutputStream()));
                obj.stdout = BufferedReader(InputStreamReader(obj.proc.getInputStream()));

                % Read initial prompt
                [~, obj.startupLines] = readToNextInputPrompt(obj);
                if isempty(obj.startupLines)
                    warning('Was expecting output from Maxima but did not get any.')
                end
            catch ME
                error('Maxima could not be started: %s', ME.message);
            end
        end

        function s = readAll(obj)
        % readAll Read all available characters from obj.stdout
        %   s = readAll(obj) reads while data is available and returns
        %   a MATLAB char vector.

            sb = javaObject('java.lang.StringBuilder');
            
            while obj.stdout.ready()
                code = obj.stdout.read();
                if code < 0
                    break;
                end
                sb.append(char(code));
            end
            s = char(sb.toString());    % convert Java String to MATLAB char
            if ~isempty(obj.transcript_file)
                fprintf(obj.transcript_fid, '%s', s);
            end
        end
    end
end
