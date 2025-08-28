---
header-includes:
    - '\usepackage[a4paper]{geometry}'
    - '\usepackage{algorithm}'
    - '\usepackage[noEnd=true]{algpseudocodex}'
documentclass: article
fontsize: 12pt
---

\begin{algorithm}[H]
    \caption{Distanza secondo Spearman}
    \label{alg:spearman_distance}
    \begin{algorithmic}[1]
		\Function{SpearmanDistance}{lista1, lista2}
			\State $d \gets 0$
			\ForAll{$gene \in$ lista1}
				\State $tmp \gets$ trova $gene$ in lista2
				\If{$\exists~tmp$}
					$d \gets d + |gene_{rank} - tmp_{rank}|$
				\EndIf
			\EndFor
			\State \textbf{return} $d$
		\EndFunction
    \end{algorithmic}
\end{algorithm}
