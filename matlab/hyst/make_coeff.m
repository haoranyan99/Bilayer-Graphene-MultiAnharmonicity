function coef = make_coeff(par)

    coef = struct();

    % --- artificial chi ---
    chi = make_chi_artificial(par);

    % --- a1, b1 ---
    a1 = @(T,dop) par.el.a1_scaler .* ...
        (par.el.a1_invV - abs(chi(T,dop)));

    b1 = @(T,dop) par.el.b1;

    % --- builders ---
    coef.build_el  = @(T,dop) ...
        build_electronic_energy_poly(struct( ...
            'a1', a1(T,dop), ...
            'b1', b1(T,dop)));

    coef.build_lat = build_lattice_potential_simple(par);

    % expose
    coef.chi = chi;
    coef.a1  = a1;
    coef.b1  = b1;
end