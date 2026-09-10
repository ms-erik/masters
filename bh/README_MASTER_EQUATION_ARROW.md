# CINCO_camadas_equacao_mestra

Esta pasta é uma nova abordagem para o ARROW eletrônico de cinco camadas.
Ela preserva o pipeline de transporte original, mas substitui a perturbação
fenomenológica imposta diretamente por uma cadeia física:

```text
m_s preparado
  -> equação mestra spin-carga NV-/NV0
  -> Delta P_NVminus(m_s)
  -> Delta V0(m_s) = e^2 n_NV DeltaP_NVminus / C_A
  -> Delta V(z,m_s)
  -> T(E), DeltaT(E), corrente de Landauer
```

## Script principal

Execute:

```matlab
run_multipair_master_equation_arrow_demo
```

Ele faz:

1. carrega os parâmetros da estrutura de cinco camadas;
2. resolve a equação mestra para os estados A e B;
3. sobrescreve `params.spin.deltaV_A_meV` e `params.spin.deltaV_B_meV`;
4. chama o mesmo transporte usado na versão anterior;
5. salva `Resultados_Multipair_MasterEquation_ARROW.mat`.

## Arquivos novos

- `nv_charge_master_matrix.m`: matriz de taxas do modelo compacto NV-/NV0.
- `nv_initial_population.m`: população inicial após preparação de spin.
- `solve_nv_charge_master.m`: resolve pulso finito ou estado estacionário.
- `update_spin_deltaV_from_master_equation.m`: converte população de carga em `DeltaV0`.
- `plot_nv_master_equation_summary.m`: gráficos diagnósticos da equação mestra.
- `run_multipair_master_equation_arrow_demo.m`: fluxo completo master equation + ARROW.
- `run_master_equation_parameter_sweep.m`: varredura rápida de densidade de NV e pulso.

## Observações importantes

Os valores de taxas incluídos em `multipair_electronic_arrow_params.m` são
valores iniciais plausíveis para simulação, não constantes universais. A
etapa experimental ou de comparação com literatura deve calibrar:

- taxas de ionização/recombinação;
- duração e potência do pulso de leitura;
- densidade areal de NV;
- capacitância efetiva por unidade de área `C_A`;
- convenção de sinal de energia potencial eletrônica.

Por padrão, a referência de potencial é o ponto médio entre os estados A e
B. Assim o dispositivo `none` continua representando a heteroestrutura sem
perturbação NV, enquanto os estados A/B aparecem como deslocamentos opostos
em torno desse baseline.
