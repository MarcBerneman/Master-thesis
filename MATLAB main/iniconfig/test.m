clear
close all

ini = IniConfig();
ini.ReadFile('example.ini');

sections = ini.GetSections();
[keys, count_keys] = ini.GetKeys(sections{1});
values = ini.GetValues(sections{1}, keys);