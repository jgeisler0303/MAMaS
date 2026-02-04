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
        maxWait = 2.0; % seconds
    end

    properties (SetAccess = private)
        startupLines
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
        function obj = MaximaInterface(maxWait)
            arguments
                maxWait (1,1) double {mustBePositive} = 2.0
            end

            obj.nextId = 1;
            obj.maxWait = maxWait;
            obj.startMaxima();
        end

        function id = sendNoWait(obj, cmd)
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

            obj.stdin.write([cmd, newline]);
            obj.stdin.flush();

            % Maxima updates the internal %o number even without output
            id = sprintf('%%o%d', obj.nextId);
            obj.nextId = obj.nextId + 1;
        end

        function result = sendAndParse(obj, cmd)
            % Send a command to Maxima and parse the response
            % Returns result string and output ID (%oN)
            if isempty(obj.proc) || ~obj.proc.isAlive()
                delete(obj);
                error('Maxima process has died. Interface and all expressions are now invalid.');
            end

            if ~endsWith(strtrim(cmd), ';')
                cmd = [cmd, ';'];
            end
            obj.stdin.write([cmd, newline]);
            obj.stdin.flush();

            [result, ~] = obj.readToNextInputPrompt();

            if isempty(result)
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
        end

        function checkAlive(obj)
            if isempty(obj.proc) || ~obj.proc.isAlive()
                error('MaximaExpression: The Maxima process has died. Expression is invalid.');
            end
        end            
        
    end

    methods (Static)
        function obj = getInstance(maxWait)
            arguments
                maxWait (1,1) double {mustBePositive} = 2.0
            end

            % Singleton accessor
            persistent instance
            if isempty(instance) || ~isvalid(instance)
                instance = MAMaS.MaximaInterface(maxWait);
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
        function [outputLines, extraLines] = readToNextInputPrompt(obj)
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

                    tok = regexp(line, '^\(%i(\d+)\)', 'once', 'tokens');
                    if ~isempty(tok)
                        lastSeenIn = str2double(tok{1});
                        % We may have to catch up to the current input id
                        if lastSeenIn > obj.nextId
                            error('Maxima input prompt IDs are not sequential. Read %%i%d but while current should be %%i%d.', lastSeenIn, obj.nextId);
                        end
                        continue;
                    end

                     [tok, id_end] = regexp(line, '^\(%o(\d+)\)', 'once', 'tokens', 'end');
                    if ~isempty(tok)
                        outputLines{end+1} = strtrim(line(id_end+1:end)); %#ok<AGROW>
                        % Output ids may be more than one at a time. Next input id must be one more than last output id.
                        obj.nextId = str2double(tok{1}) + 1;
                    else
                        extraLines{end+1} = line; %#ok<AGROW>
                    end
                end
                if lastSeenIn == obj.nextId
                    break;
                end
                pause(0.01);
            end
            if toc(tStart) >= obj.maxWait
                warning('Timeout while reading from Maxima process. Expression IDs may be out of sync.');
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
        end
    end
end
