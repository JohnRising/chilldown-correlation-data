# When is the line chilled down? (t_cd = wall reaches T_sat + 10 K)

'When is the line chilled down?' — chilldown-completion metrics for every closure and station.

  .venv/bin/python -m ofchill.chilldown     -> results/final11/chilldown_{stations.csv, summary.csv, time.md}

Definitions (same for model and data):
  t_cd     first time the wall temperature falls below T_sat + 10 K.  T_sat is the measured local saturation
           temperature for UF (station Tsat(t)) and F2 (Tsat(P_exit(t))); Tsat(P_out) for the Darr 1-g tests;
           Tsat(P_drive) for NBS (local line pressure not recorded).
  t_q      quench time = time of the maximum cooling rate (data: 5-sample running mean; model: 3-point), only
           counted when the trace drops by > 30 K within +-1.5 s of it (a real quench).
  outcome  hit (both complete within the data record), miss (data complete, model not), false (model
           completes, data not), none (neither).  The model is integrated at least to the end of the data
           record; model completion is judged within the model horizon (t_cd of a 'miss' is censored).
  whole section: last-station t_cd per test (all stations of the test must complete in the data).


## Summary per dataset and closure

| dataset     | closure                                   |   n_stations |   n_data_complete |   hit |   miss |   false |   MAE_s |   bias_s |   MAPE_pct |   within_1s |   within_2s |   within_10pct |   quench_MAE_s |   whole_n |   whole_missed |   whole_MAPE_pct |   whole_bias_pct |
|:------------|:------------------------------------------|-------------:|------------------:|------:|-------:|--------:|--------:|---------:|-----------:|------------:|------------:|---------------:|---------------:|----------:|---------------:|-----------------:|-----------------:|
| Darr 1 g    | Darr 2016 x1.19 + R5                      |            4 |                 4 |     4 |      0 |       0 |    1.81 |     1.81 |      10.82 |        0.00 |        0.50 |           0.50 |           2.33 |         2 |              0 |            11.61 |            11.61 |
| Darr 1 g    | Dittus-Boelter (sat. vapour)              |            4 |                 4 |     4 |      0 |       0 |    9.82 |    -9.82 |      58.14 |        0.00 |        0.00 |           0.00 |          15.88 |         2 |              0 |            57.17 |           -57.17 |
| Darr 1 g    | Darr 2016 film law + R5                   |            4 |                 4 |     1 |      3 |       0 |    4.38 |     4.38 |      25.57 |        0.00 |        0.00 |           0.00 |           4.70 |         2 |              2 |           nan    |           nan    |
| Darr 1 g    | Miropolskii                               |            2 |                 2 |     0 |      2 |       0 |  nan    |   nan    |     nan    |        0.00 |        0.00 |           0.00 |         nan    |         1 |              1 |           nan    |           nan    |
| Darr 1 g    | present closure                           |            4 |                 4 |     4 |      0 |       0 |    0.81 |    -0.39 |       4.84 |        0.75 |        1.00 |           0.75 |           0.98 |         2 |              0 |             2.82 |             2.22 |
| Darr 1 g    | z-independent film law + R5               |            4 |                 4 |     4 |      0 |       0 |    0.55 |    -0.40 |       3.15 |        0.75 |        1.00 |           1.00 |           0.60 |         2 |              0 |             2.45 |            -0.75 |
| F2 (0.05 g) | Darr 2016 x1.19 + R5                      |            1 |                 1 |     1 |      0 |       0 |    0.27 |     0.27 |       1.93 |        1.00 |        1.00 |           1.00 |           0.65 |         1 |              0 |             1.93 |             1.93 |
| F2 (0.05 g) | Dittus-Boelter (sat. vapour)              |            1 |                 1 |     1 |      0 |       0 |   10.22 |   -10.22 |      74.52 |        0.00 |        0.00 |           0.00 |          13.10 |         1 |              0 |            74.52 |           -74.52 |
| F2 (0.05 g) | Darr 2016 film law + R5                   |            1 |                 1 |     1 |      0 |       0 |    3.44 |     3.44 |      25.08 |        0.00 |        0.00 |           0.00 |           3.80 |         1 |              0 |            25.08 |            25.08 |
| F2 (0.05 g) | Miropolskii                               |            1 |                 1 |     0 |      1 |       0 |  nan    |   nan    |     nan    |        0.00 |        0.00 |           0.00 |         nan    |         1 |              1 |           nan    |           nan    |
| F2 (0.05 g) | present closure                           |            1 |                 1 |     1 |      0 |       0 |    0.31 |     0.31 |       2.30 |        1.00 |        1.00 |           1.00 |           0.70 |         1 |              0 |             2.30 |             2.30 |
| F2 (0.05 g) | z-independent film law + R5               |            1 |                 1 |     1 |      0 |       0 |    0.61 |    -0.61 |       4.47 |        1.00 |        1.00 |           1.00 |           0.45 |         1 |              0 |             4.47 |            -4.47 |
| NBS LH2     | Miropolskii                               |            4 |                 4 |     4 |      0 |       0 |    4.64 |    -4.64 |      10.94 |        0.25 |        0.25 |           0.50 |           3.45 |         1 |              0 |             4.39 |            -4.39 |
| NBS LH2     | present closure                           |           24 |                24 |    24 |      0 |       0 |    2.98 |    -1.99 |       6.22 |        0.25 |        0.50 |           0.88 |           4.11 |         6 |              0 |             4.96 |             0.06 |
| NBS LN2     | Darr 2016 x1.19 (contact-Tc rewet)        |            8 |                 7 |     7 |      0 |       0 |   13.86 |     9.19 |      12.56 |        0.14 |        0.14 |           0.43 |         nan    |         1 |              0 |             8.55 |             8.55 |
| NBS LN2     | Miropolskii                               |           12 |                10 |    10 |      0 |       0 |   73.40 |    73.40 |      80.69 |        0.00 |        0.00 |           0.20 |         nan    |         1 |              0 |             9.10 |             9.10 |
| NBS LN2     | present closure                           |           12 |                10 |    10 |      0 |       0 |   16.24 |     9.58 |      13.94 |        0.00 |        0.00 |           0.30 |         nan    |         1 |              0 |             8.77 |             8.77 |
| NBS LN2     | z-independent film law (contact-Tc rewet) |           12 |                10 |    10 |      0 |       0 |   16.93 |     9.77 |      14.35 |        0.00 |        0.00 |           0.30 |         nan    |         1 |              0 |             8.37 |             8.37 |
| UF          | Darr 2016 x1.19 (contact-Tc rewet)        |           20 |                 8 |     7 |      1 |       5 |    2.02 |    -0.49 |      13.75 |        0.38 |        0.50 |           0.50 |           2.76 |         4 |              1 |             8.74 |            -4.67 |
| UF          | Darr 2016 x1.19 + R5                      |           20 |                 8 |     8 |      0 |       5 |    2.30 |    -1.36 |      13.72 |        0.25 |        0.50 |           0.62 |           1.48 |         4 |              0 |             9.48 |            -6.39 |
| UF          | Dittus-Boelter (sat. vapour)              |           20 |                 8 |     8 |      0 |       8 |   12.06 |   -12.06 |      67.62 |        0.00 |        0.00 |           0.00 |          15.66 |         4 |              0 |            69.56 |           -69.56 |
| UF          | Darr 2016 (contact-Tc rewet)              |           20 |                 8 |     4 |      4 |       0 |    2.02 |     2.02 |      10.91 |        0.25 |        0.25 |           0.25 |           7.51 |         4 |              3 |             3.32 |             3.32 |
| UF          | Darr 2016 film law + R5                   |           20 |                 8 |     6 |      2 |       0 |    1.75 |     1.30 |       9.39 |        0.25 |        0.50 |           0.50 |           4.62 |         4 |              1 |            10.31 |             8.23 |
| UF          | Miropolskii                               |           18 |                 8 |     0 |      8 |       0 |  nan    |   nan    |     nan    |        0.00 |        0.00 |           0.00 |         nan    |         4 |              4 |           nan    |           nan    |
| UF          | present closure                           |           20 |                 8 |     8 |      0 |       2 |    2.09 |    -1.43 |      12.27 |        0.25 |        0.50 |           0.50 |           1.20 |         4 |              0 |             6.51 |            -6.51 |
| UF          | z-independent film law + R5               |           20 |                 8 |     8 |      0 |       2 |    2.19 |    -1.33 |      13.57 |        0.38 |        0.50 |           0.50 |           1.48 |         4 |              0 |            10.80 |           -10.80 |

## Per station

| dataset     | test   | station   | closure                                   |   t_end_data |   tcd_d |   tcd_m |   tq_d |   tq_m | outcome   |    err |   rel_err |   q_err |
|:------------|:-------|:----------|:------------------------------------------|-------------:|--------:|--------:|-------:|-------:|:----------|-------:|----------:|--------:|
| UF          | 1      | up        | present closure                           |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | dn        | present closure                           |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | up        | z-independent film law + R5               |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | dn        | z-independent film law + R5               |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | up        | Darr 2016 film law + R5                   |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | dn        | Darr 2016 film law + R5                   |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | up        | Darr 2016 x1.19 + R5                      |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | dn        | Darr 2016 x1.19 + R5                      |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | up        | Darr 2016 (contact-Tc rewet)              |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | dn        | Darr 2016 (contact-Tc rewet)              |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | up        | Darr 2016 x1.19 (contact-Tc rewet)        |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | dn        | Darr 2016 x1.19 (contact-Tc rewet)        |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | up        | Dittus-Boelter (sat. vapour)              |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | dn        | Dittus-Boelter (sat. vapour)              |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | up        | Miropolskii                               |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 1      | dn        | Miropolskii                               |        23.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | up        | present closure                           |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | dn        | present closure                           |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | up        | z-independent film law + R5               |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | dn        | z-independent film law + R5               |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | up        | Darr 2016 film law + R5                   |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | dn        | Darr 2016 film law + R5                   |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | up        | Darr 2016 x1.19 + R5                      |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | dn        | Darr 2016 x1.19 + R5                      |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | up        | Darr 2016 (contact-Tc rewet)              |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | dn        | Darr 2016 (contact-Tc rewet)              |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | up        | Darr 2016 x1.19 (contact-Tc rewet)        |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | dn        | Darr 2016 x1.19 (contact-Tc rewet)        |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | up        | Dittus-Boelter (sat. vapour)              |        20.64 |  nan    |  nan    | nan    |   1.10 | none      | nan    |    nan    |  nan    |
| UF          | 2      | dn        | Dittus-Boelter (sat. vapour)              |        20.64 |  nan    |  nan    | nan    |   1.10 | none      | nan    |    nan    |  nan    |
| UF          | 2      | up        | Miropolskii                               |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 2      | dn        | Miropolskii                               |        20.64 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 3      | up        | present closure                           |        29.36 |  nan    |  nan    |  27.45 | nan    | none      | nan    |    nan    |  nan    |
| UF          | 3      | dn        | present closure                           |        29.36 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 3      | up        | z-independent film law + R5               |        29.36 |  nan    |  nan    |  27.45 | nan    | none      | nan    |    nan    |  nan    |
| UF          | 3      | dn        | z-independent film law + R5               |        29.36 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 3      | up        | Darr 2016 film law + R5                   |        29.36 |  nan    |  nan    |  27.45 | nan    | none      | nan    |    nan    |  nan    |
| UF          | 3      | dn        | Darr 2016 film law + R5                   |        29.36 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 3      | up        | Darr 2016 x1.19 + R5                      |        29.36 |  nan    |   25.93 |  27.45 |  25.70 | false     | nan    |    nan    |   -1.75 |
| UF          | 3      | dn        | Darr 2016 x1.19 + R5                      |        29.36 |  nan    |   26.66 | nan    |  26.55 | false     | nan    |    nan    |  nan    |
| UF          | 3      | up        | Darr 2016 (contact-Tc rewet)              |        29.36 |  nan    |  nan    |  27.45 | nan    | none      | nan    |    nan    |  nan    |
| UF          | 3      | dn        | Darr 2016 (contact-Tc rewet)              |        29.36 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 3      | up        | Darr 2016 x1.19 (contact-Tc rewet)        |        29.36 |  nan    |   27.15 |  27.45 |  26.90 | false     | nan    |    nan    |   -0.55 |
| UF          | 3      | dn        | Darr 2016 x1.19 (contact-Tc rewet)        |        29.36 |  nan    |   27.89 | nan    |  27.80 | false     | nan    |    nan    |  nan    |
| UF          | 3      | up        | Dittus-Boelter (sat. vapour)              |        29.36 |  nan    |   15.99 |  27.45 |   3.30 | false     | nan    |    nan    |  -24.15 |
| UF          | 3      | dn        | Dittus-Boelter (sat. vapour)              |        29.36 |  nan    |   16.48 | nan    |   3.30 | false     | nan    |    nan    |  nan    |
| UF          | 3      | up        | Miropolskii                               |        29.36 |  nan    |  nan    |  27.45 | nan    | none      | nan    |    nan    |  nan    |
| UF          | 3      | dn        | Miropolskii                               |        29.36 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | up        | present closure                           |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | dn        | present closure                           |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | up        | z-independent film law + R5               |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | dn        | z-independent film law + R5               |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | up        | Darr 2016 film law + R5                   |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | dn        | Darr 2016 film law + R5                   |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | up        | Darr 2016 x1.19 + R5                      |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | dn        | Darr 2016 x1.19 + R5                      |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | up        | Darr 2016 (contact-Tc rewet)              |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | dn        | Darr 2016 (contact-Tc rewet)              |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | up        | Darr 2016 x1.19 (contact-Tc rewet)        |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | dn        | Darr 2016 x1.19 (contact-Tc rewet)        |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | up        | Dittus-Boelter (sat. vapour)              |        20.20 |  nan    |   12.76 | nan    |   2.10 | false     | nan    |    nan    |  nan    |
| UF          | 4      | dn        | Dittus-Boelter (sat. vapour)              |        20.20 |  nan    |   13.34 | nan    |   2.10 | false     | nan    |    nan    |  nan    |
| UF          | 4      | up        | Miropolskii                               |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 4      | dn        | Miropolskii                               |        20.20 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | up        | present closure                           |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | dn        | present closure                           |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | up        | z-independent film law + R5               |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | dn        | z-independent film law + R5               |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | up        | Darr 2016 film law + R5                   |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | dn        | Darr 2016 film law + R5                   |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | up        | Darr 2016 x1.19 + R5                      |        24.82 |  nan    |   23.38 | nan    |  23.10 | false     | nan    |    nan    |  nan    |
| UF          | 5      | dn        | Darr 2016 x1.19 + R5                      |        24.82 |  nan    |  nan    | nan    |  24.70 | none      | nan    |    nan    |  nan    |
| UF          | 5      | up        | Darr 2016 (contact-Tc rewet)              |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | dn        | Darr 2016 (contact-Tc rewet)              |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | up        | Darr 2016 x1.19 (contact-Tc rewet)        |        24.82 |  nan    |   24.44 | nan    |  24.20 | false     | nan    |    nan    |  nan    |
| UF          | 5      | dn        | Darr 2016 x1.19 (contact-Tc rewet)        |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | up        | Dittus-Boelter (sat. vapour)              |        24.82 |  nan    |   12.46 | nan    |   0.25 | false     | nan    |    nan    |  nan    |
| UF          | 5      | dn        | Dittus-Boelter (sat. vapour)              |        24.82 |  nan    |   12.85 | nan    |   0.30 | false     | nan    |    nan    |  nan    |
| UF          | 5      | up        | Miropolskii                               |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 5      | dn        | Miropolskii                               |        24.82 |  nan    |  nan    | nan    | nan    | none      | nan    |    nan    |  nan    |
| UF          | 6      | up        | present closure                           |        30.73 |  nan    |   23.62 |  23.27 |  23.15 | false     | nan    |    nan    |   -0.12 |
| UF          | 6      | dn        | present closure                           |        30.73 |  nan    |   27.35 |  29.09 |  27.25 | false     | nan    |    nan    |   -1.84 |
| UF          | 6      | up        | z-independent film law + R5               |        30.73 |  nan    |   24.68 |  23.27 |  24.60 | false     | nan    |    nan    |    1.33 |
| UF          | 6      | dn        | z-independent film law + R5               |        30.73 |  nan    |   25.85 |  29.09 |  25.45 | false     | nan    |    nan    |   -3.64 |
| UF          | 6      | up        | Darr 2016 film law + R5                   |        30.73 |  nan    |  nan    |  23.27 | nan    | none      | nan    |    nan    |  nan    |
| UF          | 6      | dn        | Darr 2016 film law + R5                   |        30.73 |  nan    |  nan    |  29.09 | nan    | none      | nan    |    nan    |  nan    |
| UF          | 6      | up        | Darr 2016 x1.19 + R5                      |        30.73 |  nan    |   26.47 |  23.27 |  26.30 | false     | nan    |    nan    |    3.03 |
| UF          | 6      | dn        | Darr 2016 x1.19 + R5                      |        30.73 |  nan    |   28.45 |  29.09 |  28.30 | false     | nan    |    nan    |   -0.79 |
| UF          | 6      | up        | Darr 2016 (contact-Tc rewet)              |        30.73 |  nan    |  nan    |  23.27 | nan    | none      | nan    |    nan    |  nan    |
| UF          | 6      | dn        | Darr 2016 (contact-Tc rewet)              |        30.73 |  nan    |  nan    |  29.09 | nan    | none      | nan    |    nan    |  nan    |
| UF          | 6      | up        | Darr 2016 x1.19 (contact-Tc rewet)        |        30.73 |  nan    |   27.84 |  23.27 |  27.65 | false     | nan    |    nan    |    4.38 |
| UF          | 6      | dn        | Darr 2016 x1.19 (contact-Tc rewet)        |        30.73 |  nan    |   30.13 |  29.09 |  30.10 | false     | nan    |    nan    |    1.01 |
| UF          | 6      | up        | Dittus-Boelter (sat. vapour)              |        30.73 |  nan    |   13.03 |  23.27 |   3.40 | false     | nan    |    nan    |  -19.87 |
| UF          | 6      | dn        | Dittus-Boelter (sat. vapour)              |        30.73 |  nan    |   13.46 |  29.09 |   3.40 | false     | nan    |    nan    |  -25.69 |
| UF          | 6      | up        | Miropolskii                               |        30.73 |  nan    |  nan    |  23.27 | nan    | no run    | nan    |    nan    |  nan    |
| UF          | 6      | dn        | Miropolskii                               |        30.73 |  nan    |  nan    |  29.09 | nan    | no run    | nan    |    nan    |  nan    |
| UF          | 7      | up        | present closure                           |        28.20 |   21.82 |   17.25 |  18.40 |  16.95 | hit       |  -4.56 |     -0.21 |   -1.45 |
| UF          | 7      | dn        | present closure                           |        28.20 |   24.19 |   19.97 |  21.13 |  19.85 | hit       |  -4.22 |     -0.17 |   -1.28 |
| UF          | 7      | up        | z-independent film law + R5               |        28.20 |   21.82 |   18.38 |  18.40 |  18.30 | hit       |  -3.43 |     -0.16 |   -0.10 |
| UF          | 7      | dn        | z-independent film law + R5               |        28.20 |   24.19 |   19.07 |  21.13 |  19.05 | hit       |  -5.12 |     -0.21 |   -2.08 |
| UF          | 7      | up        | Darr 2016 film law + R5                   |        28.20 |   21.82 |   21.21 |  18.40 |  20.90 | hit       |  -0.61 |     -0.03 |    2.50 |
| UF          | 7      | dn        | Darr 2016 film law + R5                   |        28.20 |   24.19 |   23.44 |  21.13 |  23.40 | hit       |  -0.75 |     -0.03 |    2.27 |
| UF          | 7      | up        | Darr 2016 x1.19 + R5                      |        28.20 |   21.82 |   17.29 |  18.40 |  17.10 | hit       |  -4.53 |     -0.21 |   -1.30 |
| UF          | 7      | dn        | Darr 2016 x1.19 + R5                      |        28.20 |   24.19 |   19.02 |  21.13 |  19.00 | hit       |  -5.16 |     -0.21 |   -2.13 |
| UF          | 7      | up        | Darr 2016 (contact-Tc rewet)              |        28.20 |   21.82 |   22.65 |  18.40 |  22.40 | hit       |   0.84 |      0.04 |    4.00 |
| UF          | 7      | dn        | Darr 2016 (contact-Tc rewet)              |        28.20 |   24.19 |   24.99 |  21.13 |  24.95 | hit       |   0.80 |      0.03 |    3.82 |
| UF          | 7      | up        | Darr 2016 x1.19 (contact-Tc rewet)        |        28.20 |   21.82 |   18.46 |  18.40 |  18.20 | hit       |  -3.35 |     -0.15 |   -0.20 |
| UF          | 7      | dn        | Darr 2016 x1.19 (contact-Tc rewet)        |        28.20 |   24.19 |   20.36 |  21.13 |  20.35 | hit       |  -3.83 |     -0.16 |   -0.78 |
| UF          | 7      | up        | Dittus-Boelter (sat. vapour)              |        28.20 |   21.82 |    6.46 |  18.40 |   0.85 | hit       | -15.35 |     -0.70 |  -17.55 |
| UF          | 7      | dn        | Dittus-Boelter (sat. vapour)              |        28.20 |   24.19 |    6.82 |  21.13 |   0.85 | hit       | -17.37 |     -0.72 |  -20.28 |
| UF          | 7      | up        | Miropolskii                               |        28.20 |   21.82 |  nan    |  18.40 | nan    | miss      | nan    |    nan    |  nan    |
| UF          | 7      | dn        | Miropolskii                               |        28.20 |   24.19 |  nan    |  21.13 | nan    | miss      | nan    |    nan    |  nan    |
| UF          | 8      | up        | present closure                           |        23.07 |   18.94 |   16.92 |  17.20 |  16.70 | hit       |  -2.01 |     -0.11 |   -0.50 |
| UF          | 8      | dn        | present closure                           |        23.07 |   20.81 |   19.15 |  20.67 |  19.10 | hit       |  -1.65 |     -0.08 |   -1.57 |
| UF          | 8      | up        | z-independent film law + R5               |        23.07 |   18.94 |   17.85 |  17.20 |  17.80 | hit       |  -1.09 |     -0.06 |    0.60 |
| UF          | 8      | dn        | z-independent film law + R5               |        23.07 |   20.81 |   18.34 |  20.67 |  18.35 | hit       |  -2.46 |     -0.12 |   -2.32 |
| UF          | 8      | up        | Darr 2016 film law + R5                   |        23.07 |   18.94 |   20.60 |  17.20 |  20.40 | hit       |   1.67 |      0.09 |    3.20 |
| UF          | 8      | dn        | Darr 2016 film law + R5                   |        23.07 |   20.81 |   22.38 |  20.67 |  22.30 | hit       |   1.57 |      0.08 |    1.63 |
| UF          | 8      | up        | Darr 2016 x1.19 + R5                      |        23.07 |   18.94 |   17.22 |  17.20 |  17.05 | hit       |  -1.71 |     -0.09 |   -0.15 |
| UF          | 8      | dn        | Darr 2016 x1.19 + R5                      |        23.07 |   20.81 |   18.75 |  20.67 |  18.70 | hit       |  -2.06 |     -0.10 |   -1.97 |
| UF          | 8      | up        | Darr 2016 (contact-Tc rewet)              |        23.07 |   18.94 |   21.88 |  17.20 |  21.70 | hit       |   2.94 |      0.16 |    4.50 |
| UF          | 8      | dn        | Darr 2016 (contact-Tc rewet)              |        23.07 |   20.81 |  nan    |  20.67 | nan    | miss      | nan    |    nan    |  nan    |
| UF          | 8      | up        | Darr 2016 x1.19 (contact-Tc rewet)        |        23.07 |   18.94 |   18.25 |  17.20 |  18.10 | hit       |  -0.69 |     -0.04 |    0.90 |
| UF          | 8      | dn        | Darr 2016 x1.19 (contact-Tc rewet)        |        23.07 |   20.81 |   19.92 |  20.67 |  19.90 | hit       |  -0.89 |     -0.04 |   -0.77 |
| UF          | 8      | up        | Dittus-Boelter (sat. vapour)              |        23.07 |   18.94 |    6.89 |  17.20 |   0.55 | hit       | -12.05 |     -0.64 |  -16.65 |
| UF          | 8      | dn        | Dittus-Boelter (sat. vapour)              |        23.07 |   20.81 |    7.23 |  20.67 |   0.60 | hit       | -13.58 |     -0.65 |  -20.07 |
| UF          | 8      | up        | Miropolskii                               |        23.07 |   18.94 |  nan    |  17.20 | nan    | miss      | nan    |    nan    |  nan    |
| UF          | 8      | dn        | Miropolskii                               |        23.07 |   20.81 |  nan    |  20.67 | nan    | miss      | nan    |    nan    |  nan    |
| UF          | 9      | up        | present closure                           |        21.80 |   16.78 |   15.24 |  14.07 |  15.05 | hit       |  -1.54 |     -0.09 |    0.98 |
| UF          | 9      | dn        | present closure                           |        21.80 |   17.77 |   17.67 |  17.67 |  17.60 | hit       |  -0.10 |     -0.01 |   -0.07 |
| UF          | 9      | up        | z-independent film law + R5               |        21.80 |   16.78 |   16.34 |  14.07 |  16.25 | hit       |  -0.44 |     -0.03 |    2.18 |
| UF          | 9      | dn        | z-independent film law + R5               |        21.80 |   17.77 |   16.96 |  17.67 |  16.95 | hit       |  -0.80 |     -0.05 |   -0.72 |
| UF          | 9      | up        | Darr 2016 film law + R5                   |        21.80 |   16.78 |   19.10 |  14.07 |  18.90 | hit       |   2.32 |      0.14 |    4.83 |
| UF          | 9      | dn        | Darr 2016 film law + R5                   |        21.80 |   17.77 |   21.36 |  17.67 |  21.30 | hit       |   3.60 |      0.20 |    3.63 |
| UF          | 9      | up        | Darr 2016 x1.19 + R5                      |        21.80 |   16.78 |   15.67 |  14.07 |  15.50 | hit       |  -1.11 |     -0.07 |    1.43 |
| UF          | 9      | dn        | Darr 2016 x1.19 + R5                      |        21.80 |   17.77 |   17.68 |  17.67 |  17.65 | hit       |  -0.09 |     -0.01 |   -0.02 |
| UF          | 9      | up        | Darr 2016 (contact-Tc rewet)              |        21.80 |   16.78 |   20.30 |  14.07 |  20.15 | hit       |   3.52 |      0.21 |    6.08 |
| UF          | 9      | dn        | Darr 2016 (contact-Tc rewet)              |        21.80 |   17.77 |  nan    |  17.67 |   2.40 | miss      | nan    |    nan    |  -15.27 |
| UF          | 9      | up        | Darr 2016 x1.19 (contact-Tc rewet)        |        21.80 |   16.78 |   16.73 |  14.07 |  16.55 | hit       |  -0.05 |     -0.00 |    2.48 |
| UF          | 9      | dn        | Darr 2016 x1.19 (contact-Tc rewet)        |        21.80 |   17.77 |   18.85 |  17.67 |  18.80 | hit       |   1.08 |      0.06 |    1.13 |
| UF          | 9      | up        | Dittus-Boelter (sat. vapour)              |        21.80 |   16.78 |    5.34 |  14.07 |   0.30 | hit       | -11.44 |     -0.68 |  -13.77 |
| UF          | 9      | dn        | Dittus-Boelter (sat. vapour)              |        21.80 |   17.77 |    5.63 |  17.67 |   0.30 | hit       | -12.14 |     -0.68 |  -17.37 |
| UF          | 9      | up        | Miropolskii                               |        21.80 |   16.78 |  nan    |  14.07 | nan    | miss      | nan    |    nan    |  nan    |
| UF          | 9      | dn        | Miropolskii                               |        21.80 |   17.77 |  nan    |  17.67 | nan    | miss      | nan    |    nan    |  nan    |
| UF          | 10     | up        | present closure                           |        14.27 |    8.41 |   11.04 |   8.13 |  10.95 | hit       |   2.64 |      0.31 |    2.82 |
| UF          | 10     | dn        | present closure                           |        14.27 |   12.98 |   12.96 |  12.00 |  12.95 | hit       |  -0.01 |     -0.00 |    0.95 |
| UF          | 10     | up        | z-independent film law + R5               |        14.27 |    8.41 |   11.87 |   8.13 |  11.80 | hit       |   3.46 |      0.41 |    3.67 |
| UF          | 10     | dn        | z-independent film law + R5               |        14.27 |   12.98 |   12.24 |  12.00 |  11.85 | hit       |  -0.74 |     -0.06 |   -0.15 |
| UF          | 10     | up        | Darr 2016 film law + R5                   |        14.27 |    8.41 |  nan    |   8.13 |   0.60 | miss      | nan    |    nan    |   -7.53 |
| UF          | 10     | dn        | Darr 2016 film law + R5                   |        14.27 |   12.98 |  nan    |  12.00 |   0.65 | miss      | nan    |    nan    |  -11.35 |
| UF          | 10     | up        | Darr 2016 x1.19 + R5                      |        14.27 |    8.41 |   11.38 |   8.13 |  11.25 | hit       |   2.98 |      0.35 |    3.12 |
| UF          | 10     | dn        | Darr 2016 x1.19 + R5                      |        14.27 |   12.98 |   13.78 |  12.00 |  13.75 | hit       |   0.80 |      0.06 |    1.75 |
| UF          | 10     | up        | Darr 2016 (contact-Tc rewet)              |        14.27 |    8.41 |  nan    |   8.13 |   0.60 | miss      | nan    |    nan    |   -7.53 |
| UF          | 10     | dn        | Darr 2016 (contact-Tc rewet)              |        14.27 |   12.98 |  nan    |  12.00 |   0.65 | miss      | nan    |    nan    |  -11.35 |
| UF          | 10     | up        | Darr 2016 x1.19 (contact-Tc rewet)        |        14.27 |    8.41 |   12.67 |   8.13 |  12.55 | hit       |   4.27 |      0.51 |    4.42 |
| UF          | 10     | dn        | Darr 2016 x1.19 (contact-Tc rewet)        |        14.27 |   12.98 |  nan    |  12.00 |   0.60 | miss      | nan    |    nan    |  -11.40 |
| UF          | 10     | up        | Dittus-Boelter (sat. vapour)              |        14.27 |    8.41 |    3.32 |   8.13 |   0.25 | hit       |  -5.09 |     -0.61 |   -7.88 |
| UF          | 10     | dn        | Dittus-Boelter (sat. vapour)              |        14.27 |   12.98 |    3.52 |  12.00 |   0.30 | hit       |  -9.46 |     -0.73 |  -11.70 |
| UF          | 10     | up        | Miropolskii                               |        14.27 |    8.41 |  nan    |   8.13 | nan    | miss      | nan    |    nan    |  nan    |
| UF          | 10     | dn        | Miropolskii                               |        14.27 |   12.98 |  nan    |  12.00 | nan    | miss      | nan    |    nan    |  nan    |
| F2 (0.05 g) | F2     | exit      | present closure                           |        19.70 |   13.71 |   14.03 |  13.30 |  14.00 | hit       |   0.31 |      0.02 |    0.70 |
| Darr 1 g    | 31     | up        | present closure                           |        22.85 |   17.11 |   15.29 |  16.55 |  14.95 | hit       |  -1.83 |     -0.11 |   -1.60 |
| Darr 1 g    | 31     | dn        | present closure                           |        22.85 |   18.10 |   17.99 |  17.20 |  17.95 | hit       |  -0.11 |     -0.01 |    0.75 |
| Darr 1 g    | 40     | up        | present closure                           |        20.35 |   15.47 |   15.00 |  14.90 |  14.70 | hit       |  -0.47 |     -0.03 |   -0.20 |
| Darr 1 g    | 40     | dn        | present closure                           |        20.35 |   16.78 |   17.63 |  16.20 |  17.55 | hit       |   0.85 |      0.05 |    1.35 |
| F2 (0.05 g) | F2     | exit      | z-independent film law + R5               |        19.70 |   13.71 |   13.10 |  13.30 |  12.85 | hit       |  -0.61 |     -0.04 |   -0.45 |
| Darr 1 g    | 31     | up        | z-independent film law + R5               |        22.85 |   17.11 |   15.80 |  16.55 |  15.60 | hit       |  -1.31 |     -0.08 |   -0.95 |
| Darr 1 g    | 31     | dn        | z-independent film law + R5               |        22.85 |   18.10 |   17.52 |  17.20 |  17.45 | hit       |  -0.58 |     -0.03 |    0.25 |
| Darr 1 g    | 40     | up        | z-independent film law + R5               |        20.35 |   15.47 |   15.47 |  14.90 |  15.30 | hit       |  -0.00 |     -0.00 |    0.40 |
| Darr 1 g    | 40     | dn        | z-independent film law + R5               |        20.35 |   16.78 |   17.07 |  16.20 |  17.00 | hit       |   0.29 |      0.02 |    0.80 |
| F2 (0.05 g) | F2     | exit      | Darr 2016 film law + R5                   |        19.70 |   13.71 |   17.15 |  13.30 |  17.10 | hit       |   3.44 |      0.25 |    3.80 |
| Darr 1 g    | 31     | up        | Darr 2016 film law + R5                   |        22.85 |   17.11 |   21.49 |  16.55 |  21.25 | hit       |   4.38 |      0.26 |    4.70 |
| Darr 1 g    | 31     | dn        | Darr 2016 film law + R5                   |        22.85 |   18.10 |  nan    |  17.20 | nan    | miss      | nan    |    nan    |  nan    |
| Darr 1 g    | 40     | up        | Darr 2016 film law + R5                   |        20.35 |   15.47 |  nan    |  14.90 | nan    | miss      | nan    |    nan    |  nan    |
| Darr 1 g    | 40     | dn        | Darr 2016 film law + R5                   |        20.35 |   16.78 |  nan    |  16.20 | nan    | miss      | nan    |    nan    |  nan    |
| F2 (0.05 g) | F2     | exit      | Darr 2016 x1.19 + R5                      |        19.70 |   13.71 |   13.98 |  13.30 |  13.95 | hit       |   0.27 |      0.02 |    0.65 |
| Darr 1 g    | 31     | up        | Darr 2016 x1.19 + R5                      |        22.85 |   17.11 |   18.23 |  16.55 |  18.00 | hit       |   1.12 |      0.07 |    1.45 |
| Darr 1 g    | 31     | dn        | Darr 2016 x1.19 + R5                      |        22.85 |   18.10 |   19.82 |  17.20 |  19.75 | hit       |   1.72 |      0.10 |    2.55 |
| Darr 1 g    | 40     | up        | Darr 2016 x1.19 + R5                      |        20.35 |   15.47 |   17.57 |  14.90 |  17.40 | hit       |   2.10 |      0.14 |    2.50 |
| Darr 1 g    | 40     | dn        | Darr 2016 x1.19 + R5                      |        20.35 |   16.78 |   19.08 |  16.20 |  19.00 | hit       |   2.30 |      0.14 |    2.80 |
| F2 (0.05 g) | F2     | exit      | Dittus-Boelter (sat. vapour)              |        19.70 |   13.71 |    3.49 |  13.30 |   0.20 | hit       | -10.22 |     -0.75 |  -13.10 |
| Darr 1 g    | 31     | up        | Dittus-Boelter (sat. vapour)              |        22.85 |   17.11 |    6.34 |  16.55 |   0.35 | hit       | -10.77 |     -0.63 |  -16.20 |
| Darr 1 g    | 31     | dn        | Dittus-Boelter (sat. vapour)              |        22.85 |   18.10 |    7.68 |  17.20 |   0.35 | hit       | -10.41 |     -0.58 |  -16.85 |
| Darr 1 g    | 40     | up        | Dittus-Boelter (sat. vapour)              |        20.35 |   15.47 |    6.92 |  14.90 |   0.30 | hit       |  -8.55 |     -0.55 |  -14.60 |
| Darr 1 g    | 40     | dn        | Dittus-Boelter (sat. vapour)              |        20.35 |   16.78 |    7.25 |  16.20 |   0.35 | hit       |  -9.53 |     -0.57 |  -15.85 |
| F2 (0.05 g) | F2     | exit      | Miropolskii                               |        19.70 |   13.71 |  nan    |  13.30 | nan    | miss      | nan    |    nan    |  nan    |
| Darr 1 g    | 31     | up        | Miropolskii                               |        22.85 |   17.11 |  nan    |  16.55 | nan    | no run    | nan    |    nan    |  nan    |
| Darr 1 g    | 31     | dn        | Miropolskii                               |        22.85 |   18.10 |  nan    |  17.20 | nan    | no run    | nan    |    nan    |  nan    |
| Darr 1 g    | 40     | up        | Miropolskii                               |        20.35 |   15.47 |  nan    |  14.90 | nan    | miss      | nan    |    nan    |  nan    |
| Darr 1 g    | 40     | dn        | Miropolskii                               |        20.35 |   16.78 |  nan    |  16.20 | nan    | miss      | nan    |    nan    |  nan    |
| NBS LN2     | 10     | s1        | present closure                           |       109.80 |   84.49 |   75.11 | nan    | nan    | hit       |  -9.38 |     -0.11 |  nan    |
| NBS LN2     | 10     | s2        | present closure                           |       240.00 |  178.90 |  188.29 | nan    | nan    | hit       |   9.38 |      0.05 |  nan    |
| NBS LN2     | 10     | s3        | present closure                           |       224.80 |  223.15 |  259.16 | nan    | nan    | hit       |  36.01 |      0.16 |  nan    |
| NBS LN2     | 10     | s4        | present closure                           |       249.80 |  nan    |  278.45 | nan    | nan    | none      | nan    |    nan    |  nan    |
| NBS LN2     | 10     | s1        | Miropolskii                               |       109.80 |   84.49 |  343.23 | nan    | nan    | hit       | 258.73 |      3.06 |  nan    |
| NBS LN2     | 10     | s2        | Miropolskii                               |       240.00 |  178.90 |  252.46 | nan    | nan    | hit       |  73.55 |      0.41 |  nan    |
| NBS LN2     | 10     | s3        | Miropolskii                               |       224.80 |  223.15 |  275.62 | nan    | nan    | hit       |  52.47 |      0.24 |  nan    |
| NBS LN2     | 10     | s4        | Miropolskii                               |       249.80 |  nan    |  301.85 | nan    | nan    | none      | nan    |    nan    |  nan    |
| NBS LN2     | 10     | s1        | z-independent film law (contact-Tc rewet) |       109.80 |   84.49 |   73.40 | nan    | nan    | hit       | -11.09 |     -0.13 |  nan    |
| NBS LN2     | 10     | s2        | z-independent film law (contact-Tc rewet) |       240.00 |  178.90 |  188.18 | nan    | nan    | hit       |   9.28 |      0.05 |  nan    |
| NBS LN2     | 10     | s3        | z-independent film law (contact-Tc rewet) |       224.80 |  223.15 |  259.10 | nan    | nan    | hit       |  35.95 |      0.16 |  nan    |
| NBS LN2     | 10     | s4        | z-independent film law (contact-Tc rewet) |       249.80 |  nan    |  277.97 | nan    | nan    | none      | nan    |    nan    |  nan    |
| NBS LN2     | 11     | s1        | present closure                           |       129.60 |   79.31 |   61.14 | nan    | nan    | hit       | -18.17 |     -0.23 |  nan    |
| NBS LN2     | 11     | s2        | present closure                           |       130.00 |  128.23 |  153.40 | nan    | nan    | hit       |  25.17 |      0.20 |  nan    |
| NBS LN2     | 11     | s3        | present closure                           |       174.60 |  nan    |  209.17 | nan    | nan    | none      | nan    |    nan    |  nan    |
| NBS LN2     | 11     | s4        | present closure                           |       239.80 |  199.58 |  220.93 | 193.00 | nan    | hit       |  21.36 |      0.11 |  nan    |
| NBS LN2     | 11     | s1        | Miropolskii                               |       129.60 |   79.31 |  222.01 | nan    | nan    | hit       | 142.70 |      1.80 |  nan    |
| NBS LN2     | 11     | s2        | Miropolskii                               |       130.00 |  128.23 |  196.63 | nan    | nan    | hit       |  68.39 |      0.53 |  nan    |
| NBS LN2     | 11     | s3        | Miropolskii                               |       174.60 |  nan    |  215.80 | nan    | nan    | none      | nan    |    nan    |  nan    |
| NBS LN2     | 11     | s4        | Miropolskii                               |       239.80 |  199.58 |  228.00 | 193.00 | nan    | hit       |  28.42 |      0.14 |  nan    |
| NBS LN2     | 11     | s1        | z-independent film law (contact-Tc rewet) |       129.60 |   79.31 |   60.16 | nan    | nan    | hit       | -19.15 |     -0.24 |  nan    |
| NBS LN2     | 11     | s2        | z-independent film law (contact-Tc rewet) |       130.00 |  128.23 |  155.42 | nan    | nan    | hit       |  27.19 |      0.21 |  nan    |
| NBS LN2     | 11     | s3        | z-independent film law (contact-Tc rewet) |       174.60 |  nan    |  215.20 | nan    | nan    | none      | nan    |    nan    |  nan    |
| NBS LN2     | 11     | s4        | z-independent film law (contact-Tc rewet) |       239.80 |  199.58 |  225.98 | 193.00 | nan    | hit       |  26.40 |      0.13 |  nan    |
| NBS LN2     | 11     | s1        | Darr 2016 x1.19 (contact-Tc rewet)        |       129.60 |   79.31 |   63.73 | nan    | nan    | hit       | -15.58 |     -0.20 |  nan    |
| NBS LN2     | 11     | s2        | Darr 2016 x1.19 (contact-Tc rewet)        |       130.00 |  128.23 |  154.61 | nan    | nan    | hit       |  26.38 |      0.21 |  nan    |
| NBS LN2     | 11     | s3        | Darr 2016 x1.19 (contact-Tc rewet)        |       174.60 |  nan    |  208.42 | nan    | nan    | none      | nan    |    nan    |  nan    |
| NBS LN2     | 11     | s4        | Darr 2016 x1.19 (contact-Tc rewet)        |       239.80 |  199.58 |  220.76 | 193.00 | nan    | hit       |  21.18 |      0.11 |  nan    |
| NBS LN2     | 12     | s1        | present closure                           |       129.80 |   44.54 |   38.83 | nan    | nan    | hit       |  -5.71 |     -0.13 |  nan    |
| NBS LN2     | 12     | s2        | present closure                           |        99.55 |   74.36 |   93.05 | nan    | nan    | hit       |  18.68 |      0.25 |  nan    |
| NBS LN2     | 12     | s3        | present closure                           |       129.80 |  115.66 |  123.70 | nan    | nan    | hit       |   8.04 |      0.07 |  nan    |
| NBS LN2     | 12     | s4        | present closure                           |       129.80 |  119.12 |  129.56 | 112.75 | nan    | hit       |  10.44 |      0.09 |  nan    |
| NBS LN2     | 12     | s1        | Miropolskii                               |       129.80 |   44.54 |  100.57 | nan    | nan    | hit       |  56.03 |      1.26 |  nan    |
| NBS LN2     | 12     | s2        | Miropolskii                               |        99.55 |   74.36 |  108.93 | nan    | nan    | hit       |  34.56 |      0.46 |  nan    |
| NBS LN2     | 12     | s3        | Miropolskii                               |       129.80 |  115.66 |  123.93 | nan    | nan    | hit       |   8.27 |      0.07 |  nan    |
| NBS LN2     | 12     | s4        | Miropolskii                               |       129.80 |  119.12 |  129.96 | 112.75 | nan    | hit       |  10.84 |      0.09 |  nan    |
| NBS LN2     | 12     | s1        | z-independent film law (contact-Tc rewet) |       129.80 |   44.54 |   38.99 | nan    | nan    | hit       |  -5.56 |     -0.12 |  nan    |
| NBS LN2     | 12     | s2        | z-independent film law (contact-Tc rewet) |        99.55 |   74.36 |   91.53 | nan    | nan    | hit       |  17.17 |      0.23 |  nan    |
| NBS LN2     | 12     | s3        | z-independent film law (contact-Tc rewet) |       129.80 |  115.66 |  123.21 | nan    | nan    | hit       |   7.55 |      0.07 |  nan    |
| NBS LN2     | 12     | s4        | z-independent film law (contact-Tc rewet) |       129.80 |  119.12 |  129.08 | 112.75 | nan    | hit       |   9.96 |      0.08 |  nan    |
| NBS LN2     | 12     | s1        | Darr 2016 x1.19 (contact-Tc rewet)        |       129.80 |   44.54 |   43.76 | nan    | nan    | hit       |  -0.78 |     -0.02 |  nan    |
| NBS LN2     | 12     | s2        | Darr 2016 x1.19 (contact-Tc rewet)        |        99.55 |   74.36 |   88.86 | nan    | nan    | hit       |  14.49 |      0.19 |  nan    |
| NBS LN2     | 12     | s3        | Darr 2016 x1.19 (contact-Tc rewet)        |       129.80 |  115.66 |  124.11 | nan    | nan    | hit       |   8.45 |      0.07 |  nan    |
| NBS LN2     | 12     | s4        | Darr 2016 x1.19 (contact-Tc rewet)        |       129.80 |  119.12 |  129.31 | 112.75 | nan    | hit       |  10.19 |      0.09 |  nan    |
| NBS LH2     | 2      | s1        | present closure                           |        27.94 |   26.14 |   19.51 | nan    |  17.50 | hit       |  -6.63 |     -0.25 |  nan    |
| NBS LH2     | 2      | s2        | present closure                           |        44.82 |   42.94 |   41.92 |  35.58 |  40.00 | hit       |  -1.02 |     -0.02 |    4.42 |
| NBS LH2     | 2      | s3        | present closure                           |        65.84 |   65.56 |   57.11 |  58.36 |  55.50 | hit       |  -8.45 |     -0.13 |   -2.86 |
| NBS LH2     | 2      | s4        | present closure                           |        79.96 |   65.85 |   62.83 |  64.56 |  62.00 | hit       |  -3.02 |     -0.05 |   -2.56 |
| NBS LH2     | 2      | s1        | Miropolskii                               |        27.94 |   26.14 |   19.72 | nan    |  17.50 | hit       |  -6.43 |     -0.25 |  nan    |
| NBS LH2     | 2      | s2        | Miropolskii                               |        44.82 |   42.94 |   42.09 |  35.58 |  40.00 | hit       |  -0.85 |     -0.02 |    4.42 |
| NBS LH2     | 2      | s3        | Miropolskii                               |        65.84 |   65.56 |   57.18 |  58.36 |  55.00 | hit       |  -8.39 |     -0.13 |   -3.36 |
| NBS LH2     | 2      | s4        | Miropolskii                               |        79.96 |   65.85 |   62.96 |  64.56 |  62.00 | hit       |  -2.89 |     -0.04 |   -2.56 |
| NBS LH2     | 3      | s1        | present closure                           |       150.00 |   39.44 |   39.49 | nan    |  36.00 | hit       |   0.04 |      0.00 |  nan    |
| NBS LH2     | 3      | s2        | present closure                           |       150.00 |   94.90 |   88.30 | nan    |  84.50 | hit       |  -6.60 |     -0.07 |  nan    |
| NBS LH2     | 3      | s3        | present closure                           |       150.00 |  129.72 |  121.14 | 128.25 | 117.00 | hit       |  -8.57 |     -0.07 |  -11.25 |
| NBS LH2     | 3      | s4        | present closure                           |       144.75 |  139.46 |  133.32 | 137.62 | 131.00 | hit       |  -6.14 |     -0.04 |   -6.62 |
| NBS LH2     | 4      | s1        | present closure                           |        59.84 |   21.83 |   23.11 |   9.46 |  20.50 | hit       |   1.28 |      0.06 |   11.04 |
| NBS LH2     | 4      | s2        | present closure                           |        59.84 |   49.43 |   50.30 |  40.48 |  47.50 | hit       |   0.87 |      0.02 |    7.02 |
| NBS LH2     | 4      | s3        | present closure                           |        74.80 |   71.95 |   68.53 |  68.42 |  65.50 | hit       |  -3.43 |     -0.05 |   -2.92 |
| NBS LH2     | 4      | s4        | present closure                           |        74.98 |   72.56 |   75.31 |  71.46 |  73.50 | hit       |   2.75 |      0.04 |    2.04 |
| NBS LH2     | 5      | s1        | present closure                           |        49.84 |   18.65 |   17.67 |   5.62 |  15.00 | hit       |  -0.98 |     -0.05 |    9.38 |
| NBS LH2     | 5      | s2        | present closure                           |        49.90 |   39.78 |   37.24 |  32.52 |  35.00 | hit       |  -2.53 |     -0.06 |    2.48 |
| NBS LH2     | 5      | s3        | present closure                           |        64.92 |   56.57 |   50.44 |  52.38 |  48.00 | hit       |  -6.13 |     -0.11 |   -4.38 |
| NBS LH2     | 5      | s4        | present closure                           |        64.80 |   58.67 |   55.33 |  57.54 |  54.00 | hit       |  -3.34 |     -0.06 |   -3.54 |
| NBS LH2     | 6      | s1        | present closure                           |        41.91 |   12.99 |   13.25 |   4.40 |  10.50 | hit       |   0.26 |      0.02 |    6.10 |
| NBS LH2     | 6      | s2        | present closure                           |        41.90 |   26.56 |   27.88 |  22.43 |  25.00 | hit       |   1.32 |      0.05 |    2.57 |
| NBS LH2     | 6      | s3        | present closure                           |        41.91 |   39.07 |   37.42 |  33.77 |  35.50 | hit       |  -1.65 |     -0.04 |    1.73 |
| NBS LH2     | 6      | s4        | present closure                           |        43.97 |   38.88 |   41.09 |  37.70 |  40.00 | hit       |   2.21 |      0.06 |    2.30 |
| NBS LH2     | 7      | s1        | present closure                           |        27.91 |    8.46 |    9.21 |   6.33 |   6.50 | hit       |   0.75 |      0.09 |    0.17 |
| NBS LH2     | 7      | s2        | present closure                           |        27.91 |   16.94 |   17.80 |  15.37 |  16.50 | hit       |   0.87 |      0.05 |    1.13 |
| NBS LH2     | 7      | s3        | present closure                           |        27.91 |   25.64 |   24.44 |  22.26 |  23.00 | hit       |  -1.20 |     -0.05 |    0.74 |
| NBS LH2     | 7      | s4        | present closure                           |        27.91 |   25.87 |   27.45 |  25.54 |  26.50 | hit       |   1.58 |      0.06 |    0.96 |
