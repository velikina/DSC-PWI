function OPT=get_options

% DEFAULT VALUES

% tolerance for CG
tol=1e-8;

% tolerance of reweighted iterations
tolRW=1e-3;

% max number of CG iterations
Ni=200;

% number of reweighted iterations
Nrw=1;

% Kaiser-Bessel window radius
L=2;
fov=2; % to conform to traditional Shepp-Logan phantom FOV

% overgridding factor
of=2;

% regularization parameter vector corresponds to [L1 L2] norms
TD1=[0 0];
TD2=[0 0];
TD3=[0 0];
TD4=[0 0];

%  'cart' vs. arbitary trajectory
traj = '';

% wordy output (0/1/2)
wordy = 1;

% Final structure
OPT=struct(...
    'tol',          tol,...
    'tolRW',        tolRW,...    
    'Ni',           Ni, ...
    'L',            L,  ...
    'fov',          fov, ...
    'of',           of,...
    'wordy',        wordy,...
    'Nrw',          Nrw,...
    'TD1',          TD1,...
    'TD2',          TD2,...
    'TD3',          TD3,...
    'TD4',          TD4,...
    'traj',         '');

