function cmt_tmt_modo_tm()
% Carrega um único arquivo TMT da estrutura de 5 camadas, calcula Lc e Lpi, 
% simula a transferência de potência E plota os perfis de campo.
    clc;
    close all;
    disp('===================================================================');
    disp('                         ANÁLISE CMT-NM                            ');
    disp('===================================================================');
    
   
 disp('Carregando dados da TMT de "results_TMT_modes.mat"...');
    
    try
        data = load("resultados_TMT_modes.mat");
        result = data.result;
    catch ME
        fprintf('Erro ao carregar o arquivo TMT: %s. Verifique o nome do arquivo.\n', ME.message);
        return;
    end
    
    num_modos = length(result.mode);
    if num_modos < 2
        error('ERRO: A análise CMT-NM requer 2 modos principais. Verifique o Gap e a busca TMT.');
    end
    
    % --- Extração de Dados ---
    all_betas = [result.mode.beta];
    [~, sorted_indices] = sort(real(all_betas), 'descend');
    
    idx_par = sorted_indices(1);
    idx_impar = sorted_indices(2);
    
    beta_par = all_betas(idx_par);
    beta_impar = all_betas(idx_impar);
    
    geo = result.geometry;
    % A camada 3 é o GAP. t_layers é a espessura. Convertendo para um.
    t_gap_um = geo.t_layers(3) * 1e6; 
    
    fprintf('\n--- Dados de Modos Normais ---\n');
    fprintf('  Gap Simulado: %.2f um (Camada 3)\n', t_gap_um);
    fprintf('  Modo Par: %.4e + i*%.4e 1/m\n', real(beta_par), imag(beta_par));
    fprintf('  Modo Ímpar: %.4e + i*%.4e 1/m\n', real(beta_impar), imag(beta_impar));
    % ==================================================================
    % === 4) CÁLCULO DE PARÂMETROS CMT ===
    % ==================================================================
    
    % Assumimos betas puramente reais para o caso ideal
    kappa_r = real((beta_par - beta_impar) / 2);
    beta_avg_r = real((beta_par + beta_impar) / 2);
    % Cálculo de Lc e Lpi
    Lc_m = pi / (2 * kappa_r);
    Lc_um = Lc_m * 1e6;
    Lpi_um = 2 * Lc_um;
    
    fprintf('\n--- Parâmetros CMT ---\n');
    fprintf('  Kappa (Re): %.4e 1/m\n', kappa_r);
    fprintf('  Beta Médio (Re): %.4e 1/m\n', beta_avg_r);
    fprintf('  Comprimento de Acoplamento (Lc): %.2f um\n', Lc_um);
    fprintf('  Comprimento de Batimento (Lpi): %.2f um\n', Lpi_um);
    % ==================================================================
    % === 5) SOLUÇÃO ODE45 (Caso Ideal: Beta real e Kappa real) ===
    % ==================================================================
    
    f_ode = @(z, a) ode_cmt_coupled_ideal(z, a, beta_avg_r, kappa_r);
    A0 = [1; 0]; 
    
    Z_max_um = max(10, 4 * Lc_um); 
    Z_span = [0, Z_max_um * 1e-6];
    
    disp('Resolvendo equações diferenciais acopladas (ODE45) no caso ideal...');
    [Z, A] = ode45(f_ode, Z_span, A0);
    
    % ==================================================================
    % === 6) RESULTADOS E PLOTAGEM DE POTÊNCIA (Batimento) ===
    % ==================================================================
    
    P1 = abs(A(:, 1)).^2;
    P2 = abs(A(:, 2)).^2;
    
    figure('Name', 'Troca de Potência CMT-NM (Caso Ideal)', 'Color', 'w');
    plot(Z * 1e6, P1, 'b-', 'LineWidth', 2, 'DisplayName', 'P1 (Guia 1)');
    hold on;
    plot(Z * 1e6, P2, 'r--', 'LineWidth', 2, 'DisplayName', 'P2 (Guia 2)');
    
    title(sprintf('Transferência de Potência CMT-NM: Gap=%.2f \\mu m', t_gap_um));
    xlabel('Comprimento de Interação z (\mu m)');
    ylabel('Potência Normalizada |a_j|^2');
    legend('show');
    grid on;
    ylim([0, 1.05]);
    
    % Resultados da Simulação
    fprintf('\n=== Resultados da Simulação ===\n');
    [P2_max, idx_max] = max(P2);
    Lc_sim = Z(idx_max) * 1e6;
    
    fprintf('Potência Máxima P2 alcançada: %.3f (em z=%.2f um)\n', P2_max, Lc_sim);

    % ==================================================================
    % === 7) PLOTAGEM DE PERFIS DE CAMPO (Superposição de Modos Normais) ===
    % ==================================================================
    
    % 7.1 Extração dos perfis de campo do arquivo TMT (Usando Ex_points)
    x_vec = result.mode(idx_par).x_points * 1e6; % Coordenadas x em um
    
    % Assumindo Ex_points como o campo transversal principal (como um modo TE)
    E_par   = result.mode(idx_par).Ex_points;
    E_impar = result.mode(idx_impar).Ex_points;
    
    % Usamos a parte real do campo (Ex) para plotagem, para modos sem perdas
    E_par_real   = real(E_par);
    E_impar_real = real(E_impar);
    
    % Normalização para evitar problemas de escala na plotagem
    E_par_norm   = E_par_real / max(abs(E_par_real));
    E_impar_norm = E_impar_real / max(abs(E_impar_real));
    
    % 7.2 Perfis de Potência (Superposição)
    % Superposição em z=0: Campo concentrado no Guia 1 (P1). Requer soma dos modos.
    E_z0 = E_par_norm + E_impar_norm;
    P_z0 = abs(E_z0).^2;
    P_z0 = P_z0 / max(P_z0); % Normaliza
    
    % Superposição em z=Lc: Campo concentrado no Guia 2 (P2). Requer subtração dos modos.
    E_zLc = E_par_norm - E_impar_norm;
    P_zLc = abs(E_zLc).^2;
    P_zLc = P_zLc / max(P_zLc); % Normaliza
    
    figure('Name', 'Perfis de Campo Transversal', 'Color', 'w');
    
    subplot(2,1,1);
    plot(x_vec, P_z0, 'k', 'LineWidth', 2);
    title(sprintf('Potência em z = 0 (Entrada Guia 1): Superposição E_{par} + E_{ímpar}'));
    xlabel('Posição Transversal x (\mu m)');
    ylabel('|E_x|^2 (Normalizado)');
    grid on;
    
    subplot(2,1,2);
    plot(x_vec, P_zLc, 'k', 'LineWidth', 2);
    title(sprintf('Potência em z = Lc = %.2f \\mu m (Saída Guia 2): Superposição E_{par} - E_{ímpar}', Lc_um));
    xlabel('Posição Transversal x (\mu m)');
    ylabel('|E_x|^2 (Normalizado)');
    grid on;

end

function da_dz = ode_cmt_coupled_ideal(~, a, beta_avg_r, kappa_r)
% Função para o ODE45 no caso simétrico e ideal (sem perdas)
    
    da_dz = [
        1i * beta_avg_r * a(1) + 1i * kappa_r * a(2);  % d(a1)/dz
        1i * beta_avg_r * a(2) + 1i * kappa_r * a(1)   % d(a2)/dz
    ];
end
