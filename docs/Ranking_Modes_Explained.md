# Ranking Modes Explained

| Mode | Behavior | Use Case |
| :--- | :--- | :--- |
| `M1` | Ranks strictly by Metric 1. | Single-metric evaluation. |
| `M1_M2` | Metric 1 is primary. Metric 2 can correct (swap) ranks if a lower-ranked method significantly outperforms a higher-ranked one in Metric 2. | Balancing performance vs. cost. |
| `M1_M3A` | Metric 1 is primary. Metric 2 acts strictly as a tie-breaker. | Tie-breaking without overriding primary results. |
| `M1_M2_M3` | Full hierarchy. M1 is primary. M2 corrects M1 (iterative). M3 applies two sub-logics: (1) **One-time correction** if M2 is neutral, and (2) **Iterative tie-breaking** if both M1 and M2 are neutral. | Complex multi-objective benchmarking. |
