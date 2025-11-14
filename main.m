\documentclass[lettersize,journal]{IEEEtran}
\usepackage{amsmath,amsfonts,amsthm}
\usepackage{algorithmic}
\usepackage{algorithm}
\usepackage{array}
\usepackage[caption=false,font=normalsize,labelfont=sf,textfont=sf]{subfig}
\usepackage{textcomp}
\usepackage{stfloats}
\usepackage{url}
\usepackage{verbatim}
\usepackage{graphicx}
\usepackage{cite}
\usepackage{tikz}
\usetikzlibrary{positioning,arrows,calc,arrows.meta,patterns,patterns.meta,decorations.markings,decorations.pathreplacing,positioning}
\tikzset{every node/.style={circle}, 
	strike through/.append style={
		decoration={markings, mark=at position 0.5 with {
				\draw[-] ++ (-2pt,-2pt) -- (2pt,2pt);}
		},postaction={decorate}}
}
\usepackage{mathrsfs}
\usepackage{mathtools}
\usepackage{dutchcal}
\usepackage{titlesec}
\titlespacing*{\section}
{0pt}{1.25ex}{0.75ex}
\titlespacing*{\subsection}
{0pt}{1ex}{0.75ex}

\usepackage{hyperref}
\usepackage{cleveref}
\usepackage{hhline}
\usepackage{relsize}

\hyphenation{}

\newcommand{\Longupdownarrow}{\Big\Updownarrow}
\newcommand{\G}{\mathcal{G}}
\newcommand{\HH}{H}
\newcommand{\TT}{T}
\newcommand{\PP}{\mathcal{P}}
\newcommand{\CC}{\mathcal{C}}
\newcommand{\C}{\mathscr{C}}
\DeclareMathOperator*{\EE}{\textit{E}}
\DeclareMathOperator*{\VV}{\textit{V}}
\newcommand{\Break}{\State \textbf{break} }

\DeclareMathOperator*{\E}{W}

\newtheorem{theorem}{Theorem}
\newtheorem{lemma}{Lemma} 
\newtheorem{problem}{Problem} 
\newtheorem{corollary}{Corollary} 
\newtheorem{assumption}{Assumption}
\newtheorem{definition}{Definition} 
\newtheorem{proposition}{Proposition} 
\newtheorem{remark}{Remark}

\def\BibTeX{{\rm B\kern-.05em{\sc i\kern-.025em b}\kern-.08em
    T\kern-.1667em\lower.7ex\hbox{E}\kern-.125emX}}

\setlength{\textfloatsep}{0.2cm}

\begin{document}

\title{Coupling strength allocation using cycle basis for synchronization of dynamical networks}

\author{Aandrew Baggio Sahaya Arokiadoss and G. Arunkumar
        % <-this % stops a space
\thanks{}% <-this % stops a space
\thanks{}}

% The paper headers
\markboth{}%
{}

\IEEEpubid{}
% Remember, if you use this you must call \IEEEpubidadjcol in the second
% column for its text to clear the IEEEpubid mark.

\maketitle

\begin{abstract}
This work proposes a novel method for synchronizing a class of diffusively coupled dynamical networks. Our approach originates from identifying a key flaw in a previous spectral graph theory–based method, which claimed that solving certain linear inequalities would yield coupling strengths ensuring synchronization. We show that this claim fails for a dynamical network if its connectivity digraph belongs to a class of digraphs and has more than four vertices. To address this, we develop a purely graph-theoretic framework using the cycle basis of each strongly connected component and an undirected star graph to derive valid synchronization conditions, removing reliance on generic inequality solvers. Under the same assumptions as the earlier work, we formally prove that the derived coupling strengths guarantee synchronization. Moreover, for a digraph with $n$ vertices, our method reduces computational effort from evaluating $\binom{n}{2}$ undirected paths to only $n-1$ directed paths, leveraging strong connectivity to achieve a runtime complexity of $\mathcal{O}(n^{3})$.
\end{abstract}

\begin{IEEEkeywords}
Coupling strength allocation, Cycle basis, coupled dynamical systems, Spectral graph theory, Synchronization.
\end{IEEEkeywords}

\label{sec:introduction}
\IEEEPARstart{S}{ynchronization} in dynamical networks is a fundamental phenomenon studied across disciplines such as neuroscience~\cite{breakspear2010generative,popovych2014control}, physics~\cite{levis2017synchronization}, biology~\cite{karakaya2022effective,gonze2005spontaneous}, and engineering~\cite{lu2002synchronization,sadaoui2011predictive}. It arises when multiple identical or nearly identical systems interact to exhibit coherent behavior. A key challenge lies in understanding how network topology and system dynamics jointly determine collective behavior. Two main frameworks address synchronization: the Master Stability Function (MSF)~\cite{pecora1998master} and Lyapunov-based methods~\cite{wu1995synchronization}. The MSF requires computing Lyapunov exponents—particularly difficult for chaotic systems—while Lyapunov methods yield matrix inequalities (MIs) defining coupling strength conditions for synchronization~\cite[Theorem 4.4]{wu2007synchronization}. In time-varying networks, every change in connectivity demands solving new MIs, increasing computational cost. When topology varies much faster than dynamics, averaging approximations apply, but these fail for slowly varying networks. To address this,~\cite{liu2015synchronization} proposed the Generalized Connection Graph Method (GCGM), which used spectral graph theory to simplify MIs into linear inequalities based on computing one undirected path between each vertex pair~\cite[Theorem~4]{liu2015synchronization}. However, GCGM has conceptual and computational limitations. It assumes pre-specified vertex imbalances that depend on yet-unknown coupling strengths, rendering the method incomplete. Moreover, it fails for digraphs containing two strongly connected components (one a directed cycle, the other a single source vertex), contradicting its claimed generality. We propose a method that resolves these issues and efficiently exploits network structure to ensure synchronization. Our main contributions are:
\begin{itemize}
\item \textbf{Identification of Algorithmic Flaw:} We show that the GCGM algorithm~\cite{liu2015synchronization} fails to produce positive coupling strengths for a particular class of digraphs.
\item \textbf{Guaranteed Construction of Coupling Strengths:} We develop a corrected approach that ensures at least one feasible set of positive coupling strengths for all digraphs.
\item \textbf{Topology-Driven Efficiency:} By leveraging strongly connected components (SCCs), we reformulate coupling strength allocation as a graph traversal problem, significantly reducing computational effort.
\end{itemize}

\section{Preliminaries}

\subsection{Graphs and Digraphs}
A \textit{digraph} is an ordered pair $(\VV, \EE)$, where $\VV(\G) = \{v_1, \dots, v_n\}$ is its set of vertices and $\EE(\G) = \{e_1, \dots, e_m\} \subseteq \VV \times \VV$ is its set of arcs of $\G$. An \textit{undirected graph} (or graph) is also an ordered pair $(\VV, \EE)$, where $\VV$ is defined the same and $\EE$ is a multiset of unordered vertex pairs.  
Throughout this work, we consider \textit{simple digraphs}, i.e., digraphs without self-loops or parallel arcs. Graphs may include parallel edges, but a \textit{simple graph} contains neither.

\begin{table}[!t]
\caption{Functions, terms and definitions}
\centering
\renewcommand{\arraystretch}{1.3}
\begin{tabular}{|p{2.5cm}|p{5cm}|}
\hline
\textbf{Term / Function} & \textbf{Definition / Description} \\
\hline
$e_k = (v_i, v_j)$ & Arc from tail ($v_i$) to head ($v_j$) \\
\hline
$e_k = \{v_i, v_j\}$ & Edge between $v_i$ and $v_j$ \\
\hline
End vertices & Vertices forming an arc (edge) \\
\hline
Incoming (Outgoing) arc & Arc with the vertex as head (tail)\\
\hline
Indegree (Outdegree) & Number of incident incoming (outgoing arcs) \\
\hline
Degree & Number of incident edges\\
\hline
Adjacent arcs & Arcs sharing at least one end vertex \\
\hline
Source vertex & Vertex with indegree $= 0$ \\
\hline
Parallel arcs (edges) & Share identical end vertices \\
\hline
Self-loops & Arcs (edges) connecting a vertex to itself \\
\hline
Underlying graph & Graph formed by replacing the arcs between each vertex pair in a digraph with an edge \\
\hline
Traversal & Ordered vertex sequence with that each consecutive pair corresponding to the end vertices of an arc in the digraph\\
\hline
$\VV(\G)$ & Vertex set of $\G$\\
\hline
$\EE(\G)$ & Arc (or edge) set of $\G$\\
\hline
$|\G|$ & Number of arcs in $\G$\\
\hline
$||\G||$ & Number of vertices in $\G$\\
\hline
$\HH$ & Maps an arc to its head\\
\hline
$\TT$ & Maps an arc to its tail\\
\hline
\end{tabular}
\label{tab:Graph definitions}
\end{table}

A \textit{weighted digraph} is a triplet $(\VV, \EE, \E)$, where $\E = [w_1, \dots, w_m]^\top \in \mathbb{R}^m$ is the vector associated with each arc weight.

\begin{definition}[Vertex Imbalance]
For vertex $v_i$ in a weighted digraph, the vertex imbalance $\text{d}_i$ is the difference between the sum of outgoing and incoming arc weights :
\[\text{d}_i = \qquad
\smashoperator{\sum_{\left\{\substack{s \in \EE(\G): \HH(e_s) = v_i}\right\}}} w_s\qquad
- \qquad\smashoperator{\sum_{\left\{\substack{t \in \EE(\G): \TT(e_t) = v_i}\right\}}} w_t\]
\end{definition}
\subsection{Undirected Subgraphs}
A \textit{subgraph} $\mathcal{H}$ of $\G$ is a graph such that $\VV(\mathcal{H}) \subseteq \VV(\G)$ and $\EE(\mathcal{H}) \subseteq \EE(\G)$, and is denoted as $\mathcal{H} \subseteq \G$.  
A \textit{path} $\PP$ is a graph with $\VV(\PP) = \{v_{k_0}, \dots, v_{k_\ell}\}$ and $\EE(\PP) = \{e_{k_1}, \dots, e_{k_\ell}\}$, where $e_{k_i} = \{v_{k_{i-1}}, v_{k_i}\}$. Two vertices in a graph are said to be \textit{connected} if there exists a path subgraph that includes both vertices. A graph is \textit{connected} if every pair of its distinct vertices is connected.  
A \textit{cycle} is a connected graph in which every vertex has degree two. A \textit{tree} is a connected graph with no cycles, and a \textit{spanning tree} is a tree subgraph that includes all the vertices of the graph.  
A \textit{complete graph}, denoted by $\mathbb{K}_n$, is one in which every vertex is adjacent to every other vertex. A \textit{star graph} is a graph with one central vertex adjacent to all other vertices, while the remaining vertices are adjacent only to the central vertex.
\subsection{Directed Subgraphs}
A \textit{directed path} $\PP$ satisfies $\VV(\PP) = \{v_{k_0}, \dots, v_{k_\ell}\}$ and $\EE(\PP) = \{e_{k_1}, \dots, e_{k_\ell}\}$ with $e_{k_i} = (v_{k_{i-1}}, v_{k_i})$.  
Its tail and head are $\TT(\PP) = \TT(e_{k_1})$ and $\HH(\PP) = \HH(e_{k_\ell})$.  
A digraph is connected if its underlying graph is connected. A \textit{directed cycle} is a digraph with at most one arc between any of its vertex pairs and its underlying graph is a cycle (arcs need not align).  
A \textit{directed spanning tree} is a digraph that has a single source vertex with a unique directed path from it to every other vertex. The \textit{length} of a directed path or cycle equals its number of arcs. A digraph is \textit{strongly connected} if a directed path exists between every ordered distinct vertex pair.
\subsection{Matrices and Vectors}
For $\text{X} \in \mathbb{R}^n$, $\text{X} > 0$ (resp. $\ge 0$) indicates all entries are positive (resp. nonnegative). The zero and all-one vectors are $\mathbf{0}_n$ and $\mathbf{1}_n$ respectively. For a weighted digraph $\G$ with $n$ vertices and $m$ arcs, the incidence ($\text{Q}$) and Laplacian ($\text{L}$) matrices are :
\begin{align*}
[\text{Q}]_{ij} &=
\begin{cases}
\phantom{-}1, & \text{ if } v_i = \TT(e_j),\\
-1, & \text{ if } v_i = \HH(e_j),\\
\phantom{-}0, & \text{otherwise,}
\end{cases} \\
[\text{L}]_{ij} &=
\begin{cases}
\;\quad\smashoperator{\sum_{\left\{\substack{e_k\in \EE(\G):\\e_k = (j,i)}\right\}}}w_k, &\text{ if } i=j,\\
\phantom{\quad \sum}-w_k, &\text{ if } i \neq j \text{ and } e_k = (i,j) \in \EE(\G),\\
\phantom{-\qquad \sum}0, & \text{otherwise.}
\end{cases}
\end{align*}

$\text{L}_{\G}$ denotes the Laplacian of $\G$. In unweighted cases, all $w_{k}=1\quad \forall 1\leq k \leq m$. The vertex imbalance vector D satisfies Q$\E$=D, where the $i^{\text{th}}$ entry of D is d$_i$. Each directed cycle $\mathcal{H}$ in $\G$ has a \textit{signed incidence vector} $\text{X} = [x_1, \dots, x_m]^\top \in \mathbb{R}^{||\G||}$:
\begin{equation*}
x_k =
\begin{cases}
\phantom{-}1, & e_k \in \EE(\mathcal{H}) \text{ traversed tail to head},\\
-1, & e_k \in \EE(\mathcal{H}) \text{ traversed head to tail},\\
\phantom{-}0, & e_k \notin \EE(\mathcal{H}).
\end{cases}
\end{equation*}

Traversal orientation determines sign convention, though it is irrelevant for cycles from directed ear decompositions.  
The \textit{cycle space} $\C$ of $\G$ is the nullspace of Q : $\C = \{\text{X} \in \mathbb{R}^m \mid \text{QX} = \mathbf{0}_n\}$. A \textit{cycle basis} is a set of linearly independent signed incidence vectors spanning $\C$ and corresponding to distinct directed cycles.

\section{Problem Formulation}
We consider a network of $n$ diffusively coupled dynamical systems following \cite{liu2015synchronization} :
\begin{equation}
\dot{\boldsymbol{z}}_{i} = f(\boldsymbol{z}_{i}) + \sum w_{k}(t)\,\text{P}(\boldsymbol{z}_{j} - \boldsymbol{z}_{i})
\label{eq:Dynamical Network Model}
\end{equation}
The summation is over all the other dynamical systems that influence the system $i$. Here, $f$ is the system dynamics, $\boldsymbol{z}_i \in \mathbb{R}^d$ the state, $w_k(t)$ is the coupling strength of $e_{k}$ at time $t$ and $\text{P}$ is a $\{0,1\}$ diagonal matrix indicating which state components are coupled. Assuming the connectivity remains fixed for some time period, we omit the time dependence of $w_k$. The couplings and their strengths form a weighted digraph $\G$, termed the \textit{connectivity digraph}, where the arc weights represent coupling strengths. The synchronization manifold of~\eqref{eq:Dynamical Network Model} is defined as the set $\left\{\, \mathbf{z} \;\middle|\; \mathbf{z} = \mathbf{1}_n \otimes \boldsymbol{z},\ \boldsymbol{z} \in \mathbb{R}^d \,\right\},$ that is, the network is said to be \emph{synchronized} if its state asymptotically converges to the trajectories within this set. Global stability of this manifold is ensured under the following assumptions :
\begin{enumerate}
    \item There exists an $a>0$ such that $\forall\, \boldsymbol{x},\boldsymbol{y} \in \mathbb{R}^d$,
    \begin{equation}
    (\boldsymbol{x} - \boldsymbol{y})^{\top}[f(\boldsymbol{x}) - f(\boldsymbol{y}) - a\text{P}(\boldsymbol{x} - \boldsymbol{y})] \leq 0
    \label{eq:Synchronization Assumption}
    \end{equation}
    \item The connectivity digraph contains a directed spanning tree.
\end{enumerate}

The objective is to find an arc weight vector $\E$ satisfying 
\[
\lim_{t \to \infty} \|\boldsymbol{z}_i(t) - \boldsymbol{z}(t)\|_2 = 0,\quad \forall\, i=1,\dots,n.
\]
A Lyapunov-based approach leads to the matrix inequality \cite[Theorem~4.4]{wu2007synchronization}:
\begin{equation}
(\text{U} \otimes \text{V}) \big( \text{L}_{\G} \otimes (-\text{P}) - \text{I}_{n} \otimes \text{Y} \big) \preceq 0,
\label{eq:Synchronization Matrix Inequality}
\end{equation}
which can be simplified to be \cite{liu2013coupling}:
\begin{equation}
\text{L}_{\G} \text{L}_{\G_0} - a \text{L}_{\G_0} \succeq 0.
\label{eq:Graph inequality}
\end{equation}
For $\text{L}_{\G_0}=\text{L}_{\text{K}_n}$, the synchronization condition becomes
\begin{equation}
w_k^{s} > \frac{a}{n} b_k, \quad \forall\, k=1,\dots,m
\label{eq:Synchronization Condition}
\end{equation}
Given $e_k=(i,j)\in\EE(\G)$, $w_k^{s}=w_k+w_l$ if $e_l=(j,i)\in\EE(\G)$, and $w_k^{s}=w_k$ otherwise. Also, 
$b_k =\quad \smashoperator{\sum_{\left\{\substack{\PP \in \mathcal{S}:\\e_k \in \EE(\PP)}\right\}}} \Omega(\PP)$,  
with $\mathcal{S}$ denoting a set of undirected paths, one for each unordered distinct vertex pair of $\G$ and having the same as its end vertices in the corresponding underlying graph. $\Omega(\PP)$ is the weight associated with each such path, defined as :
\[
\Omega(\PP) =
\begin{cases}
|\PP|\chi\!\left(1+\dfrac{\text{d}_{\HH(\PP)}+\text{d}_{\TT(\PP)}}{2a}\right), & \text{ if }e_k \notin \EE(\G)\\
1+\dfrac{\text{d}_{\HH(e_k)}+\text{d}_{\TT(e_k)}}{2a},&  \text{ else }
\end{cases}
\]
where $\chi(z)=\max\{z,0\}$. Having outlined the method in \cite{liu2015synchronization}, we next identify its key issues.

\subsection{Error Due to Arbitrary Vertex Imbalances}
The GCGM algorithm assigns coupling strengths by processing one SCC at a time \cite[Sec.~IV]{liu2015synchronization}. For each SCC, a subgraph is formed, and a new vertex with outgoing arcs to vertices receiving inter-SCC arcs is added. The arc weights computed for these arcs are redistributed among existing arcs with the same head. SCCs without incoming arcs (root SCCs) are processed directly, as illustrated in \ref{fig:processed SCC}.  
This process results from the fact that we can replace multiple synchronized systems influencing one unsynchronized system by a single equivalent synchronized system as underscored in \cite{liu2015synchronization}.  
\begin{remark}
\label{thm:SCC with 1 or 2 components}
Each subgraph processed by GCGM is either:
\begin{enumerate}
\item strongly connected, or
\item composed of two SCCs, one being a single source vertex.
\end{enumerate}
\end{remark}
This processing reduces complexity by limiting path computations to smaller subgraphs. 
Consider a digraph with two SCCs: one forming a directed cycle and the other a source vertex connected to a single vertex of the cycle. 
Such a digraph arises when a directed cycle SCC has only one vertex receiving inter-SCC arcs. 
Let $w_{1}$ be the vertex imbalance of the source vertex, $-w_{1}$ be of its adjacent vertex, and $0$ for all others. 
This is achieved by assigning edge weight $w_{1}$ to the arc incident on the source vertex and $w_{2}$ to all remaining arcs. 
Using the shortest path between any two vertices and \eqref{eq:Synchronization Condition}, the synchronization inequality becomes
\begin{equation}
w_{1} >
\begin{cases}
\dfrac{2a}{n}\!\left(\text{P}^{o}_{\text{sum}} + (\text{P}^{o}_{\text{sum}} - 1)\dfrac{w_{1}}{2a}\right), & \text{ if } n~\text{ is odd},\\[6pt]
\dfrac{2a}{n}\!\left(\text{P}^{e}_{\text{sum}} + (\text{P}^{e}_{\text{sum}} - 1)\dfrac{w_{1}}{2a}\right), & \text{ if } n~\text{ is even}.
\end{cases}
\end{equation}
Here, $\text{P}^{o}_{\text{sum}} = \dfrac{n^{2} + 2n - 3}{4}$ and $\text{P}^{e}_{\text{sum}} = \dfrac{n^{2} + 2n - 4}{4}$. 
For $n > 4$, the inequality yields negative solutions, rendering positive arc weights infeasible. 
This arises from a polynomial path-sum numerator and a linear denominator. 
Feasibility can be restored by assigning vertex imbalances whose path-sums grow linearly (or slower) with $n$, a critical issue overlooked in \cite{liu2015synchronization}.

\subsection{Growth in Path Computations}
For a graph with $n$ vertices, $\frac{n(n-1)}{2}$ paths must be computed, which rapidly increases computational cost.
\subsection{Reformulated Problem}
From \cite{liu2015synchronization}:
\begin{remark}
For any path $\PP$ whose end vertices have imbalance $\leq -a$, $\chi(\Omega(\PP)) = 0$
\label{thm:Negative Vertex Imbalance}
\end{remark}
Such paths do not affect \eqref{eq:Synchronization Condition} and can be ignored. Thus, ensuring most vertex imbalances $\leq -a$ minimizes path computations. By the handshaking lemma, at most $n - 1$ vertices can have negative imbalance with positive arc weights.

\begin{lemma}[Handshaking Lemma for Weighted Digraphs]
\begin{equation*}
\sum_{v_{i}\in \VV(\G)} 
\left(
\;\;\;\;\;\;\smashoperator{\sum_{\{e_{k}\mid v_{i}=\TT(e_{k})\}}}w_{k}
\quad
-
\quad\smashoperator{\sum_{\{e_{\ell}\mid v_{i}=\HH(e_{\ell})\}}} w_{\ell}\;\;
\right)=0
\iff 
\smashoperator{\sum_{v_{i}\in \VV(\G)}} \text{d}_{i} = 0
\end{equation*}
\label{lem:handshake}
\end{lemma}

\begin{definition}[Negative Imbalance Arc Weight Vector]
A \textit{negative imbalance arc weight vector}, $\text{X}^{-}$, is a positive arc weight vector that makes $n-1$ vertices of an $n$-vertex digraph have negative imbalances.
\end{definition}

Since $\text{QX}^{-}=\text{D}$ is linear, scaling X$^{-}$ ensures all imbalances $<-a$. We thus restate the problem :
\begin{problem}
Given $a > 0$ and a weighted digraph $\G = (\VV, \EE, \E)$ that is either strongly connected or consists of a source vertex and a non-trivial SCC, find $\E > 0$ such that:
\begin{enumerate}
\item The synchronization condition \eqref{eq:Synchronization Condition} holds, and
\item All but one vertex in the non-trivial SCC have imbalance $\leq -a$.
\end{enumerate}
\end{problem}


\section{Coupling Strength Allocation}
We divide the task of determining arc weights in the connectivity digraph into intra-SCC weight allocation and inter-SCC weight allocation. For intra-SCC arcs i.e, arcs with end vertices within the same SCC, we prove that the strong connectedness allows the construction of a negative imbalance arc-weight vector using only \(n - 1\) directed paths. A feasible weight vector is then obtained as a linear combination of this imbalance vector and a vector from the SCC’s cycle space, yielding positive entries that satisfy~\eqref{eq:Modified Synchronization Condition}.  
For inter-SCC arcs, we show that the Laplacian of an $n$-vertex star graph ensures compliance with~\eqref{eq:Graph inequality}.

\subsection{Intra-SCC Arc Weight Assignment}
We first show that any strongly connected weighted digraph admits a negative imbalance arc-weight vector.

\begin{proposition}
\label{thm:maximum negative vertex imbalance vector}
Let $\PP$ be a directed path in a weighted digraph $\G$ with arc-weight vector $\E = [\,w_1, \dots, w_m\,]^{\top}$ and vertex imbalance vector D $= [\,\text{d}_1, \dots, \text{d}_n\,]^{\top}$ satisfying:
\begin{enumerate}
    \item  $w_j > w_i > 0 \;\;\forall$ $e_i, e_j \in \EE(\PP)$ with $\HH(e_j) = \TT(e_i)$ and
    \item $w_i = 0$ for all $e_i \notin \EE(\PP)$.
\end{enumerate}
Then:
\begin{enumerate}
    \item $\text{d}_i > 0$ if $v_i = \TT(\PP)$,
    \item $\text{d}_i < 0$ for $v_i \in \VV(\PP) \setminus \{\TT(\PP)\}$,
    \item $\text{d}_i = 0$ for $v_i \notin \VV(\PP)$.
\end{enumerate}
\end{proposition}

\begin{proof}
Let $\PP$ be a directed path of length $\ell$ with arcs $e_1 = (v_0, v_1), \dots, e_\ell = (v_{\ell-1}, v_\ell)$, where $v_0 = \TT(\PP)$ and $v_\ell = \HH(\PP)$.  
Assign $w_1 < \cdots < w_\ell$ on $\PP$ and $0$ elsewhere, giving
\[
\text{D} = [\,w_1,\, w_2 - w_1,\, \dots,\, -w_\ell,\, 0, \dots, 0\,]^{\top},
\]
which satisfies the stated properties.
\end{proof}

\begin{corollary}
\label{thm:negative imbalance arc weight vector for strongly connected digraphs}
Every strongly connected digraph admits a negative imbalance arc-weight vector.
\end{corollary}

\begin{proof}
Let $\G$ be strongly connected with incidence matrix Q and an arbitrary vertex $v_r$. Since $\G$ has a directed spanning tree rooted at $v_r$, each $v_j \neq v_r$ has a directed path from $v_r$.  
By Proposition~\ref{thm:maximum negative vertex imbalance vector}, each such path defines a vector X$^{-}_j$ that results in d$_j < 0$, d$_r > 0$, and d$_k \le 0$  $\forall\,k \neq r,j$. Summing over all $j \neq r$ results in a vector (X$^{-} = \sum_{j \neq r}$X$^{-}_j$), whose imbalance vector has one positive and the rest negative entries, proving the claim.
\end{proof}

To simplify computation, we set $w_i = w_j - 1 > 0$ if $\HH(e_j) = \TT(e_i)$.  
To ensure any path not starting at $v_r$ has negative weight i.e., $\Omega(\PP) < 0$, we use $\E^{-} = a\,\text{X}^{-}$ as the arc weight vector, thus changing the synchronization condition in \eqref{eq:Synchronization Condition} to be \begin{equation}
w_{k}^{s} \ge \frac{2a}{n}
\left(1 + \smashoperator{\sum_{\PP \in \mathcal{S}}} \|\PP\|\right)
\left(\smashoperator{\sum_{\PP \in \mathcal{S}}} \|\PP\|\right).
\label{eq:Modified Synchronization Condition}
\end{equation}
Since $\E^{-}$ alone may not satisfy \eqref{eq:Modified Synchronization Condition}, we add a vector from the cycle space $\mathscr{C}$ that preserves vertex imbalances i.e., Q(X$^{0} + \text{W}^-) = \text{Q}\text{W}^-, \,\, \forall\, \text{X}^{0} \in \mathscr{C}$. We construct X$^{0}$ as the sum of the nonnegative directed ear basis vectors \cite[Section~2.1]{loebl2001some} and scale it by a factor of
\[
\Delta w = \frac{2a}{n}
\left(1 + \sum_{\PP \in \mathcal{S}} \|\PP\|\right)
\left(\sum_{\PP \in \mathcal{S}} \|\PP\|\right),
\]
to get the following arc-weight vector
\[
\text{W}_{\text{sync}} = \text{W}^- + \Delta w\,\text{X}^{0} = \text{W}^- + \text{W}^0,
\]
which satisfies the synchronization condition while minimizing path computations.

\subsection{Inter-SCC Arc Weight Assignment}
After intra-SCC weight assignment, we turn towards the arcs connecting different SCCs i.e., intra-SCC arcs. For a given SCC, its \textit{upstream} SCCs are those supplying incoming arcs. After the processing as in Remark~\ref{thm:SCC with 1 or 2 components}, the resulting digraph contains two SCCs, and the following theorem helps assogn the inter-SCC weight assignment ensuring~\eqref{eq:Graph inequality}.
\begin{theorem}
\label{thm:Arc weights for digraphs with trivial SCCs}
Let $\text{L}_{\G}$ be the Laplacian of a weighted digraph $\G$ comprising one non-trivial SCC and one trivial SCC (a single source vertex). For any $a>0$, there exists an arc-weight vector such that
\begin{equation}
\frac{1}{2}\left(\text{L}_{\G_0}\text{L}_{\G} + \text{L}_{\G}^\top \text{L}_{\G_0}\right) \succeq a\,\text{L}_{\G_0},
\label{eq:symmetric-form of Graph inequality}
\end{equation}
where $\text{L}_{\G_0}$ denotes the Laplacian of an undirected star graph satisfying $\VV(\G)=\VV(\G_0)$, with its center vertex sharing the label of the source vertex of $\G$.
\end{theorem}
\begin{proof}
Let L$_{\G}$ be the weighted Laplacian matrix of a digraph of the form mentioned in Remark~\ref{thm:SCC with 1 or 2 components} having the source vertex labeled as $v_{0}$. Let L$_{\G_{0}}$ be the Laplacian matrix of an undirected star $\G_{0}$ with the same vertex set as $\G$. Let $v_0$ be the label of its central vertex. One arbitrary adjacent vertex of $v_{0}$ in $\G$ is designated the label $v_{1}$ to act as the vertex with positive vertex imbalance. By Corollary~\ref{thm:negative imbalance arc weight vector for strongly connected digraphs}, the nontrivial SCC admits a negative imbalance vector X$^{-}$ with positive imbalance at $v_1$.  
Let $\text{L}_{\text{non-triv}}$ be its weighted Laplacian matrix, and $\text{D} = [\,\text{d}_1 \cdots \text{d}_{n-1}\,]^{\top}$ its corresponding imbalance vector. Denote the outgoing arcs of $v_0$ as $\{e_{0,1}, \dots, e_{0,m_1}\}$, and define $\E_{\text{triv}} = [\,w_{0,1}\,\cdots\,w_{0,m_1}\,\boldsymbol{0}_{n-m_1-1}^{\top}]^{\top}$.  
Then, the left hand side of \ref{eq:symmetric-form of Graph inequality} can be split as follows :
\[
\begin{bmatrix}
0 & \mathbf{0}_{n-1}^\top \\[2pt]
\mathbf{0}_{n-1} &  \text{L}^{\text{sym}}_{\text{non-triv}}+\dfrac{1}{2}\text{diag}(\E)
\end{bmatrix}
+
\begin{bmatrix}
\text{d}_0 & \dfrac{1}{2}\text{D}^\top - \E_{\text{triv}}^\top \\[2pt]
\dfrac{1}{2}\text{D} - \E_{\text{triv}} & \text{diag}(\E_{\text{triv}} - \dfrac{1}{2}\text{D})
\end{bmatrix}
\]
The term L$^{\text{sym}}_{\text{non-triv}}=\dfrac{1}{2}(\text{L}_{\text{non-triv}} + \text{L}_{\text{non-triv}}^\top)$ makes the first matrix to be positive semi-definite and thus can be ignored. The second matrix share the zero non-zero pattern of a star Laplacian, and the condition further simplfies to :
%\[
%%\text{d}_0 & \dfrac{1}{2}\text{D}^\top - \E_{\text{triv}}^\top \\[2pt]
%\dfrac{1}{2}\text{D} - \E_{\text{triv}} & \text{diag}(\E_{\text{triv}} - \dfrac{1}{2}\text{D})
%\end{bmatrix}
%\succeq
%a
%\begin{bmatrix}
%n-1 & -\mathbf{1}^\top \\[2pt]
%-\mathbf{1} & \text{I}_{n}
%\end{bmatrix}
%\]
\[
w_{0,1} \ge \dfrac{\text{d}_1}{2} + a, \quad w_{0,i} \ge 0 \;\; \forall\, 1<i\le m_1
\]
Since $w_{0,i}$ are free, the inequalities always have a feasible solution.
\end{proof}
\subsection{Example}

We consider a network of Lorentz oscillators with the sequence of connectivity digraphs $\G_{k_{1}}$, $\G_{k_{2}}$, and $\G_{k_{3}}$ corresponding to the time intervals $[0,1)$, $[1,2)$, and $[2,3)$, respectively, as shown in Fig.~\ref{fig:timevarying_digraphs}. The arc weights for $\G_{k_{1}}$ are computed using our proposed method. The digraph $\G_{k_{1}}$ consists of two SCCs: $\G_{k_{1}}^{1}$ (solid blue) and $\G_{k_{1}}^{2}$ (dotted red). The upstream SCC $\G_{k_{1}}^{1}$ is condensed into vertex $v_{0}$ (Fig.~\ref{fig:processed SCC}), resulting in the reduced digraph $\tilde{\G}_{k_{1}}^{2}$ (Fig.~\ref{fig:processed SCC}).

\begin{figure}[!t]
\subfloat[]{%
\begin{tikzpicture}[>=stealth,scale=0.25, 
    redvertex/.style={circle, draw, fill=red!30, pattern=dots, pattern color=black!50!red, minimum size=14pt, inner sep=0pt},
    bluevertex/.style={circle, draw, fill=blue!30, minimum size=14pt, inner sep=0pt},
    greenvertex/.style={circle, draw, fill=green!30, pattern=north west lines, pattern color=black!40!green, minimum size=14pt, inner sep=0pt},
    arrow/.style={thick}]
    \foreach \i/\angle/\style in {
        1/90/redvertex,2/30/redvertex,3/-30/redvertex,4/-90/bluevertex,5/-150/bluevertex,6/150/bluevertex}
        \node[style=\style] (g1v\i) at (\angle:4cm) {\i};
    \node[above=of g1v1,yshift=-1cm] {$\G_{k_{1}}$};
    \draw[->,arrow] (g1v1) -- (g1v2);
    \draw[->,arrow] (g1v2) -- (g1v3);
    \draw[->,arrow] (g1v3) -- (g1v1);
    \draw[->,arrow] (g1v4) -- (g1v1);
    \draw[->,arrow] (g1v4) -- (g1v2);
    \draw[->,arrow] (g1v5) -- (g1v1);
    \draw[->,arrow] (g1v4) -- (g1v5);
    \draw[->,arrow] (g1v5) -- (g1v6);
    \draw[->,arrow] (g1v6) -- (g1v4);
\end{tikzpicture}
\label{fig:0 to 1}}
\hfil
\subfloat[]{%
\begin{tikzpicture}[>=stealth,scale=0.25, 
    redvertex/.style={circle, draw, fill=red!30, pattern=dots, pattern color=black!50!red, minimum size=14pt, inner sep=0pt},
    bluevertex/.style={circle, draw, fill=blue!30, minimum size=14pt, inner sep=0pt},
    greenvertex/.style={circle, draw, fill=green!30, pattern=north west lines, pattern color=black!40!green, minimum size=14pt, inner sep=0pt},
    arrow/.style={thick}]
    \foreach \i/\angle/\style in {
        1/90/redvertex,2/30/redvertex,3/-30/redvertex,4/-90/redvertex,5/-150/bluevertex,6/150/bluevertex}
        \node[style=\style] (g2v\i) at (\angle:4cm) {\i};
    \node[above=of g2v1,yshift=-1cm] {$\G_{k_{2}}$};
    \draw[->,arrow] (g2v1) -- (g2v2);
    \draw[->,arrow] (g2v2) -- (g2v3);
    \draw[->,arrow] (g2v3) -- (g2v4);
    \draw[->,arrow] (g2v4) -- (g2v1);
    \draw[->,arrow] (g2v5) -- (g2v6);
    \draw[->,arrow] (g2v6) -- (g2v5);
    \draw[->,arrow] (g2v4) -- (g2v5);
    \draw[->,arrow] (g2v3) -- (g2v6);
\end{tikzpicture}
\label{fig:1 to 2}}
\hfil
\subfloat[]{%
\begin{tikzpicture}[>=stealth,scale=0.25, 
    redvertex/.style={circle, draw, fill=red!30, pattern=dots, pattern color=black!50!red, minimum size=14pt, inner sep=0pt},
    bluevertex/.style={circle, draw, fill=blue!30, minimum size=14pt, inner sep=0pt},
    greenvertex/.style={circle, draw, fill=green!30, pattern=north west lines, pattern color=black!40!green, minimum size=14pt, inner sep=0pt},
    arrow/.style={thick}]
    \foreach \i/\angle/\style in {
        1/90/redvertex,2/30/redvertex,3/-30/bluevertex,4/-90/bluevertex,5/-150/greenvertex,6/150/greenvertex}
        \node[style=\style] (g3v\i) at (\angle:4cm) {\i};
    \node[above=of g3v1,yshift=-1cm] {$\G_{k_{3}}$};
    \draw[->,arrow] (g3v1) -- (g3v2);
    \draw[->,arrow] (g3v2) -- (g3v1);
    \draw[->,arrow] (g3v2) -- (g3v3);
    \draw[->,arrow] (g3v3) -- (g3v4);
    \draw[->,arrow] (g3v4) -- (g3v3);
    \draw[->,arrow] (g3v5) -- (g3v6);
    \draw[->,arrow] (g3v6) -- (g3v5);
    \draw[->,arrow] (g3v6) -- (g3v1);
\end{tikzpicture}
\label{fig:2 to 3}}
\\
\subfloat[]{%
\begin{tikzpicture}[>=stealth,scale=0.25, 
    redvertex/.style={circle, draw, fill=red!30, pattern=dots, pattern color=black!50!red, minimum size=14pt, inner sep=0pt},
    bluevertex/.style={circle, draw, fill=blue!30, minimum size=14pt, inner sep=0pt},
    greenvertex/.style={circle, draw, fill=green!30, pattern=north west lines, pattern color=black!40!green, minimum size=14pt, inner sep=0pt},
    arrow/.style={thick}]
    \foreach \i/\angle/\style in {
        1/90/redvertex,2/30/redvertex,3/-30/redvertex,0/-150/bluevertex}
        \node[style=\style] (g1v\i) at (\angle:4cm) {\i};
    \draw[->,arrow] (g1v1) -- (g1v2);
    \draw[->,arrow] (g1v2) -- (g1v3);
    \draw[->,arrow] (g1v3) -- (g1v1);
    \draw[->,arrow] (g1v5) -- (g1v2);
    \draw[->,arrow] (g1v5) -- (g1v1);
\end{tikzpicture}
\label{fig:processed SCC}}
\caption{Connectivity evolution of the dynamical network. Each color/pattern indicates a distinct SCC. Fig.~\ref{fig:processed SCC} shows the root SCC of $\G_{k_{1}}$ condensed into a single vertex.}
\label{fig:timevarying_digraphs}
\end{figure}

We compute the negative imbalance arc weight vector and directed ear basis for $\G_{k_{1}}^{1}$ (Table~\ref{tab:G1}), yielding $\Delta w = 8a$ and $\E^0 = a[\,8\,8\,8\,]^\top$. Hence, the synchronization arc weight vector for $\G_{1}$ is $\E_{\text{sync}} = a[\,11\,9\,8\,]^\top$.

\begin{table}[H]
\caption{}
\centering
\begin{tabular}{|c||c|}
\hline
\textbf{Directed Path}& \textbf{Negative Imbalance Arc Weight Vector}\\\hline
$5\rightarrow6$&$\left[\;1\;0\;0\;\right]^{\top}$\\\hline
$5\rightarrow6\rightarrow4$&$\left[\;2\;1\;0\;\right]^{\top}$\\ \hhline{|=|=|}
\textbf{Directed cycle} & \textbf{Signed Incidence Vector}\\\hline
$5\rightarrow6\rightarrow4\rightarrow5$ & $\left[\;1\;1\;1\;\right]^{\top}$\\\hline
\end{tabular}
\label{tab:G1}
\end{table}

\begin{table}[H]
\caption{}
\centering
\begin{tabular}{|c||c|}
\hline
\textbf{Directed Path}& \textbf{Negative Imbalance Arc Weight Vector}\\\hline
$1\rightarrow2$&$\left[\;1\;0\;0\;\right]^{\top}$\\\hline
$1\rightarrow2\rightarrow3$&$\left[\;2\;1\;0\;\right]^{\top}$\\\hhline{|=|=|}
\textbf{Directed cycle} & \textbf{Signed Incidence Vector}\\\hline
$1\rightarrow2\rightarrow3\rightarrow1$ & $\left[\;1\;1\;1\;\right]^{\top}$\\\hline
\end{tabular}
\label{tab:G2}
\end{table}

For $\G_{k_{1}}^{2}$, $\E^{-} = a[\,3\,1\,0\,]^\top$ (starting vertex $2$) and $X^{0} = [\,1\,1\,1\,]^{\top}$.  
To ensure a nonzero arc weight for $(3,1)$, we add the directed ear basis (Table~\ref{tab:G2}) to $\E^-$, applying zero-padding for arcs from vertex $0$.  
The final synchronization vector for $\tilde{\G}_{k_{1}}^{2}$ is $\E_{\text{sync}} = [\,3a+1\;\;a+1\;\;1\,]^\top$. Without the weights of the inter-SCC arcs, this vector gives a vertex imbalance vector D$ = a[\,3\, -2\, -1\,]^\top$. For the inter-SCC arc incident on $v_{1}$ i.e., ($0$,$1$), the arc weight needs to be $\frac{3a}{2}+a=2.5a$ and $a$ for the arc ($0$,$2$) following Theorem~\ref{thm:Arc weights for digraphs with trivial SCCs}. The weight $2.5a$ for $(0,1)$ is equally distributed between $(4,1)$ and $(5,1)$. The arc weights for the remaining digraphs are computed analogously.
\section{Simulation and Results}
We considered the Lorentz attractor system described by :\begin{equation}
\begin{aligned}
  \dot{x} &= \sigma (y - x), &
  \dot{y} &= x (r - z) - y, &
  \dot{z} &= x y - b z,&
\end{aligned}
\end{equation}
with the following parameters $\sigma = 10$, $r = 28$ and $b = 8/3$ (chaotic behavior) and computed the system parameter $a$ using the formula $a = \dfrac{b(b+1)(r+\sigma)^{2}}{16(b-1)} -\sigma$ as directed in \cite{liu2013coupling}. The initial conditions for each system were assigned randomly.
\begin{figure}[h]
\centering
    \includegraphics[scale=0.35]{Synchronization_image.eps}
    \caption{A network of 6 Lorentz oscillators, with connectivity changes as depicted in \ref{fig:timevarying_digraphs}, synchronizes via coupling strength allocation after each connectivity change.}
    \label{fig:Error Convergence}
\end{figure}
At each time step, the Euclidean norm of the difference between the $i^{\text{th}}$ state and all other states was plotted. The sum consistently decreased to zero as seen in \ref{fig:Error Convergence} indicating network synchronization irrespective of the connectivity changes.
\subsection{Complexity}
For a subgraph with $n$ vertices and $m$ arcs formed from an SCC, it involves computing $n-1$ directed paths along with determining a circuit basis. This task has the same computational complexity as finding a cycle basis using a spanning tree, as discussed in \cite{ryan1981comparison}, which is : $\mathcal{O}(n\gamma)$, where $\gamma = m - n + 1$. In the worst case scenario, where $m = \frac{n(n-1)}{2}$, the time complexity becomes $\mathcal{O}(n^3)$.
\section{Conclusion}
\label{sec:conclusion}
This study identified a critical flaw in an existing synchronization framework for diffusively coupled dynamical networks. The earlier method fails to yield valid (positive) coupling strengths for digraphs with two SCCs with one being a single source vertex and the other a directed cycle with more than four vertices—when vertex imbalances are assigned arbitrarily. To address this, we developed a graph-theoretic synchronization method that guarantees positive couplings without solving linear or matrix inequalities. The method leverages SCC structure to allocate intra and inter-component couplings systematically. By replacing exhaustive path computations with only $n-1$ directed paths, it achieves $O(n^3)$ complexity.
\section*{Acknowledgment}
The authors acknowledge the use of language assistance tools, including ChatGPT, Gemini, DeepSeek and Writefull for paraphrasing suggestions and grammar refinement during manuscript preparation. These tools were utilized solely for linguistic enhancement and not for the generation of any technical or scientific content.
 
\bibliographystyle{ieeetr}
\bibliography{References}

\end{document}


