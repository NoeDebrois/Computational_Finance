function F = CHAR_FCT(u, parameters, flag)
    % flag=0 --> Characteristic function for backward problem (with conjugate)
    % flag=1 --> Characteristic function for forward problem (no conjugate)

    % Default behaviour applies complex conjugate:
    if nargin==2
        flag=0;
    end
    
    dt = parameters.dt;
    
    % Characteristic fct for the desired model under RN-measure w/o
    % rf-rate/potential dividends:
    switch parameters.distr
    case 1 % Normal
        F = exp(dt * CHAR_EXP_BS(u,parameters));
    case 2 % Merton
        F = exp(dt * CHAR_EXP_MERTON(u,parameters));
    case 3 % Kou
        F = exp(dt * CHAR_EXP_KOU(u,parameters));
    case 4 % NIG
        F = exp(dt * CHAR_EXP_NIG(u,parameters));
    case 5 % VG 
        F = exp(dt * CHAR_EXP_VG(u,parameters));
    %case 6 % ExtVG
    %    F = exp(dt * CHAR_EXP_EXTVG(u,parameters));
    %case 7 % ExtNIG
    %    F = exp(dt * CHAR_EXP_EXTNIG(u,parameters));
    end

    % Correction to take into account the rf-rate, & potential dividends:
    F = F .* exp((parameters.rf - parameters.q) .* 1i .* u * dt);

    % Correction to apply complex conjugate if flag==0:
    if flag==0
        F = conj(F);
    end