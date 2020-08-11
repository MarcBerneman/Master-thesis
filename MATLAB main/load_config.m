function [keys,values] = load_config(section)
    section = char(section);
    ini = IniConfig();
    ini.ReadFile('configuration.ini');
    
    [keys, count_keys] = ini.GetKeys(section);
    values = ini.GetValues(section, keys);
end

