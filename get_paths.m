% Created on 3 Mar 2015
% Primary author: Kiefer Forseth
%                 kjforseth@gmail.com
%                 602.531.6430
%
% Summary
%   Find the valid directory structure on any system.
% 
% Output
%   paths: structure containing accessible directories
% 
% Example
%   paths = get_paths;
% 
% Change log
%   24 Mar 15: updated sector finding to find active Sector[0124]-[\d]

function paths = get_paths

% identify user & dropbox location
LocalHost = java.net.InetAddress.getLocalHost;
Name    = char(LocalHost.getHostName);
Address = char(LocalHost.getHostAddress);
switch Name
    case {'Arams-MacBook-Pro'}
        cloud   = '/Users/kiefer/Dropbox/tandon/std/';
        temp    = '~/Desktop/edata/std/';
        volumes = '~/Desktop/edata/Volumes/';
        code    = '~/Desktop/edata/X/';
        root    = '~/Desktop/edata/';
    case {'Arams-MacBook.local'}
        cloud   = '/Users/kiefer/Dropbox/tandon/std/';
        temp    = '~/Desktop/edata/std/';
        volumes = '~/Desktop/edata/Volumes/';
        code    = '~/Desktop/edata/X/';
        root    = '~/Desktop/edata/';
    case {'KJF-macbook.local' 'msinfected177.med.uth.tmc.edu'}
        cloud   = '/Users/kiefer/Dropbox/tandon/std/';
        temp    = '~/Desktop/tandon_tmp/';
        volumes = '/Users/kiefer/Dropbox/Volumes/';
    case 'NSURG-MBP-TANLAB1.attlocal.net'
        cloud   = '/Users/kiefer/Dropbox/tandon/std/';
        temp    = '~/Desktop/tandon_tmp/';
        volumes = '/Users/kiefer/Dropbox/Volumes/';
    case 'KJF-ASUS'
        cloud   = 'D:\Dropbox\tandon\std\';
        temp    = 'D:\tandon_tmp\';
        volumes = 'V:\';
    case 'Cyberpower'
        cloud   = 'D:\Dropbox\tandon\std\';
        temp    = 'D:\tandon_tmp\';
        volumes = 'D:\Dropbox\Volumes\';
%         volumes = 'V:\';
    case 'NSURG-TAN-MPR10.local'
        cloud   = '/Users/kieferforseth/Dropbox/tandon/std/';
        temp    = '~/Desktop/tandon_tmp/';
        volumes = '/Users/kieferforseth/Dropbox/Volumes/';
%         volumes = '/Volumes/';
    case 'NSURG-TAN-MPRO8.local' % Kamin old trashcan
        cloud   = '/Users/kiefer/Dropbox/tandon/std/';
        temp    = '~/Desktop/tandon_tmp/';
        volumes = '/Users/kiefer/Dropbox/Volumes/';
    case 'NSURG-TAN-IMAC2.local' % Kamin left
        cloud   = '/Volumes/Sector4/Users/KJF/dropbox_clone/tandon/std/';
        temp    = '~/Desktop/tandon_tmp/';
        volumes = '/Volumes/';
    case 'Mac-Pro.local' % Matt
        cloud   = '/Volumes/Sector4/Users/KJF/dropbox_clone/tandon/std/';
        temp    = '~/Desktop/tandon_tmp/';
        volumes = '/Volumes/';
    otherwise
        fprintf('Name: %s\nAddress: %s\n',Name,Address)
        error('Unrecognized user, update get_paths() for this computer')
end

% define operating system
if     ispc,   paths.os = 'pc';
elseif ismac,  paths.os = 'mac';
elseif isunix, paths.os = 'unix';
else error('Unable to determine operating system.')
end

% define relative paths for std folder
paths.std.root = temp;
paths.std.data = [paths.std.root 'data' filesep];
paths.std.docs = [paths.std.root 'docs' filesep];
paths.std.fcns = [paths.std.root 'fcns' filesep];
paths.std.figs = [paths.std.root 'figs' filesep];
paths.std.guis = [paths.std.root 'guis' filesep];
paths.std.libs = [paths.std.root 'libs' filesep];
addpath(genpath(paths.std.root))

% define relative paths for sectors
paths.vol.S0 = [];
paths.vol.S1 = [];
paths.vol.S2 = [];
paths.vol.S4 = [];

home = cd;
dirList = dir(volumes);
for d = 1:length(dirList)
    if length(dirList(d).name) < 7, continue, end
    
    paths_cur = [volumes dirList(d).name filesep];
    [~,perm,~] = fileattrib(paths_cur);
    switch dirList(d).name(1:7)
        case 'Sector0'
            if isempty(paths.vol.S0) && perm.UserWrite, paths.vol.S0 = paths_cur; end
        case 'Sector1'
            if isempty(paths.vol.S1) && perm.UserWrite, paths.vol.S1 = paths_cur; end
        case 'Sector2'
            if isempty(paths.vol.S2) && perm.UserWrite, paths.vol.S2 = paths_cur; end
       case 'Sector4'
            if isempty(paths.vol.S4) && perm.UserWrite, paths.vol.S4 = paths_cur; end
    end
end
cd(home)

% define tmp location
paths.tmp = temp;
% define code location
paths.code = code;
paths.root = root;
if ~exist(paths.tmp,'dir'), mkdir(paths.tmp), end