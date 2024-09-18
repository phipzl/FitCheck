function config = load_config(filename)
    % Read configuration file and return a struct with fields
    config = struct();
    fid = fopen(filename, 'r');
    while ~feof(fid)
        line = fgetl(fid);
        if isempty(line) || line(1) == '#'
            continue;  % Skip comments and empty lines
        end
        tokens = regexp(line, '(\w+)_PATH="([^"]+)"', 'tokens');
        if ~isempty(tokens)
            key = tokens{1}{1};
            value = tokens{1}{2};
            config.([key '_PATH']) = value;
        end
    end
    fclose(fid);
end

