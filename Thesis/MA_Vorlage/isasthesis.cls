
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Identification
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{isasthesis}[2012/03/24 The ISAS Thesis Class]


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%            Class Options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{ifthen}

\newboolean{isEnglish}
\setboolean{isEnglish}{false}
\DeclareOption{english}{\setboolean{isEnglish}{true}}
\DeclareRobustCommand{\ifenglish}[2]{\ifthenelse{\boolean{isEnglish}}{#1}{#2}}

\newboolean{hasTables}
\setboolean{hasTables}{false}
\DeclareOption{tables}{\setboolean{hasTables}{true}}
\DeclareRobustCommand{\iftables}[2]{\ifthenelse{\boolean{hasTables}}{#1}{#2}}

\newboolean{hasAlgorithms}
\setboolean{hasAlgorithms}{false}
\DeclareOption{algorithms}{\setboolean{hasAlgorithms}{true}}
\DeclareRobustCommand{\ifalgorithms}[2]{\ifthenelse{\boolean{hasAlgorithms}}{#1}{#2}}

\newboolean{hasReviews}
\setboolean{hasReviews}{false}
\DeclareOption{reviews}{\setboolean{hasReviews}{true}}
\DeclareRobustCommand{\ifreviews}[2]{\ifthenelse{\boolean{hasReviews}}{#1}{#2}}

\newboolean{hasExamples}
\setboolean{hasExamples}{false}
\DeclareOption{examples}{\setboolean{hasExamples}{true}}
\DeclareRobustCommand{\ifexamples}[2]{\ifthenelse{\boolean{hasExamples}}{#1}{#2}}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Thesis preferences
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Thesis Type
\def\thesistype#1{\gdef\@thesistype{#1}}
\global\let\@thesistype\@empty

% 1st Supervisor
\def\firstsupervisor#1{\gdef\@firstsupervisor{#1}}
\global\let\@firstsupervisor\@empty

% 2nd Supervisor
\def\secondsupervisor#1{\gdef\@secondsupervisor{#1}}
\global\let\@secondsupervisor\@empty

% Discussant
\def\discussant#1{\gdef\@discussant{#1}}
\global\let\@discussant\@empty


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Main declarations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ProcessOptions\relax

\RequirePackage{ifpdf}

\ifpdf
	\LoadClass[pdftex,twoside,openright,numbers=noenddot]{scrreprt}
\else
	\LoadClass[twoside,openright,numbers=noenddot]{scrreprt}
\fi

\RequirePackage{isastextmacros}
\RequirePackage{isasmathmacros}

\RequirePackage[chapter]{algorithm}
\RequirePackage{algorithmic}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Distinguish between PDF and DVI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifpdf
	\RequirePackage[pdftex]{graphicx}
	\RequirePackage[pdftex]{color}
	\RequirePackage[pdftex,
					pdfstartview=FitH,
					pdfpagemode={UseOutlines},
					bookmarks=true,
					bookmarksnumbered=true,
					hypertexnames=false,
					colorlinks=false,
					linkcolor=red,
					citecolor=green]{hyperref}
\else
	\RequirePackage[dvips]{graphicx,color}
	\RequirePackage[hypertex=true,
					plainpages=false,
					pdfpagelabels=true]{hyperref}
\fi

\RequirePackage[all]{hypcap} 	% Link images (instead of captions) to ref
\RequirePackage{caption} 		% Fancy captions


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           German language
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifenglish{}{
	\RequirePackage[ngerman]{babel}
	\RequirePackage{ngerman}
	\RequirePackage{bibgerm}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Various Layout Settings
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{abstract}

\renewcommand{\abstractname}{Zusammenfassung}
\renewcommand{\abstractnamefont}{\normalfont\bfseries}
\renewcommand{\abstracttextfont}{\normalfont}


\RequirePackage{fancyhdr}	% Farmatting of headers

\captionsetup{font=small,labelfont=bf}

\setlength{\parindent}{0em}
\setlength{\parskip}{1.5ex plus0.5ex minus0.2ex}
\renewcommand{\baselinestretch}{1.2}

\setlength{\fboxrule}{1pt}
\setlength{\fboxsep}{3mm}

\renewcommand{\normalsize}{\fontsize{12pt}{14pt}\selectfont}

% Specify paage margins
\newlength{\pagemarginleft}
\setlength{\pagemarginleft}{2.5cm}
\newlength{\pagemarginright}
\setlength{\pagemarginright}{3.5cm}

\RequirePackage[a4paper,
			    left=\pagemarginleft,
			    right=\pagemarginright,
			    top=2cm,
			    bottom=1.5cm,
			    includeheadfoot]{geometry}
\RequirePackage{changepage}

\setlength{\headheight}{14pt}
\pagestyle{fancy}
\renewcommand{\chaptermark}[1]{\markboth{\chaptername\ \thechapter.\ #1}{}} 
\fancyhead[LE]{\sffamily \leftmark}     			% Chapter in the left on even pages
\fancyhead[RO]{\sffamily \nouppercase{\rightmark}}	% Section in the right on odd pages
\fancyhead[RE]{}									% No text  in the right on even pages
\fancyhead[LO]{}									% No text in the left on odd pages
\fancyfoot[OR,EL]{\sffamily \thepage}
\fancyfoot[OC,EC]{\sffamily}
\renewcommand{\headrulewidth}{0.2pt}
\renewcommand{\footrulewidth}{0.2pt}

\fancypagestyle{plain}{
	\fancyhead{}
	\renewcommand{\headrulewidth}{0pt}
	\fancyfoot[RE]{}
}

\newfont{\HUGEoli}{cmssbx10 scaled4302}
\newcommand{\thesisfont}[1]{\HUGEoli #1}
\newcommand{\thesistitlefont}[1]{\huge\rmfamily\upshape\mdseries #1}
\newcommand{\namefont}[1]{\LARGE\rmfamily\upshape\bfseries #1}
\newcommand{\datefont}[1]{\Large\rmfamily\upshape\mdseries #1}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Title page
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\maketitle}{
	\thispagestyle{empty}
	% Center the titlepage (remove margins)
	\begin{adjustwidth*}{-\pagemarginleft}{-\pagemarginright}
		\begin{center}
			\begin{tabular}{>{\centering\arraybackslash}p{15cm}}
				\hline\\\includegraphics[width=15cm]{isasheader.pdf}\vspace{0.3cm}\\\hline
			\end{tabular}
			
			\vspace{2.0cm}
			
			\thesisfont{\@thesistype}
			
			\vspace{1.5cm}
			
			\thesistitlefont{\@title}
			\vfill
			
			\namefont{\@author}
			\vfill
			
			\datefont{\@date}
			\vfill
			
			\begin{tabular*}{15cm}{ll}
				\hline
				& \\
				Referent:  & Prof. Dr.-Ing. Uwe D. Hanebeck \\
				\ifthenelse{\equal{\@discussant}{}}{}{			% Show only if discussant is present
					Korreferent:  & \@discussant \\
				}
				& \\
				Betreuer:  & \@firstsupervisor \\
				\ifthenelse{\equal{\@secondsupervisor}{}}{}{	% Show only if 2nd supervisor is present
					& \@secondsupervisor \\
				}
				& \\
				\hline
			\end{tabular*}
		\end{center}
	\end{adjustwidth*}
	
	\cleardoublepage
}


\newcommand{\maketoc}{
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%%           Declaration of Academic Honesty
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	\chapter*{Eidesstattliche Erkl\"arung}
		
		\thispagestyle{empty}
		\vspace*{\fill}
		\begin{par}
			\noindent
			Hiermit erkl\"are ich, die vorliegende \@thesistype ~selbstst\"andig angefertigt zu haben. Die verwendeten Quellen sind im Text gekennzeichnet und im Literaturverzeichnis aufgef\"uhrt.
			
			\vspace{4\baselineskip}
			
			Karlsruhe, \@date
			\hspace{1cm}
			\parbox[t]{6cm}{
				\begin{center}
				\hrulefill\\ {\small \@author}
				\end{center}
			}
		\end{par}

	\cleardoublepage
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Table of Contents
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	\pagenumbering{Roman}\setcounter{page}{1}
	
	\pagestyle{fancy}
	
	\fancyhead[LE,RO]{\sffamily
		\ifenglish{
			Table of Contents
		}{
			Inhaltsverzeichnis
		}
	}
	
	\ifpdf
		\ifenglish{
			\pdfbookmark[1]{Table of Contents}{toc}
		}{
			\pdfbookmark[1]{Inhaltsverzeichnis}{toc}
		}
	\fi
	
	\tableofcontents
	
	\cleardoublepage
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           List of Figures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	\fancyhead[LE,RO]{\sffamily
		\ifenglish{
			List of Figures
		}{
			Abbildungsverzeichnis
		}
	}
	
	\phantomsection
	\ifenglish{
		\addcontentsline{toc}{chapter}{List of Figures}
	}{
		\addcontentsline{toc}{chapter}{Abbildungsverzeichnis}
	}
	
	\listoffigures 
	
	\cleardoublepage

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           List of Tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	\iftables{
		\fancyhead[LE,RO]{\sffamily
			\ifenglish{
				List of Tables
			}{
				Tabellenverzeichnis
			}
		}
	
		\phantomsection
		\ifenglish{
			\addcontentsline{toc}{chapter}{List of Tables}
		}{
			\addcontentsline{toc}{chapter}{Tabellenverzeichnis}
		}
	
		\listoftables
	
		\cleardoublepage
	}
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           List of Algorithms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	\ifalgorithms{
		\fancyhead[LE,RO]{\sffamily
			\ifenglish{
				List of Algorithms
			}{
				Algorithmenverzeichnis
			}
		}
		
		\phantomsection
		\ifenglish{
			\addcontentsline{toc}{chapter}{List of Algorithms}
		}{
			\addcontentsline{toc}{chapter}{Algorithmenverzeichnis}
		}
		
		\listofalgorithms
		
		\cleardoublepage
	}
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           List of Examples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	\ifexamples{
		\fancyhead[LE,RO]{\sffamily
			\ifenglish{
				List of Examples
			}{
				Beispielverzeichnis
			}
		}
		
		\phantomsection
		\ifenglish{
			\addcontentsline{toc}{chapter}{List of Examples}
		}{
			\addcontentsline{toc}{chapter}{Beispielverzeichnis}
		}
		
		\chapter*{
			\ifenglish{
				List of Examples
			}{
				Beispielverzeichnis
			}
		}
		
		{\sffamily \listofexamples}
		
		\cleardoublepage
	}
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           List of Review Material
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

	\ifreviews{
		\fancyhead[LE,RO]{\sffamily
			\ifenglish{
				List of Review Material
			}{
				Wiederholungsmaterialverzeichnis
			}
		}
		
		\phantomsection
		\ifenglish{
			\addcontentsline{toc}{chapter}{List of Review Material}
		}{
			\addcontentsline{toc}{chapter}{Wiederholungsmaterialverzeichnis}
		}
		
		\chapter*{
			\ifenglish{
				List of Review Material
			}{
				Wiederholungsmaterialverzeichnis
			}
		}
		
		{\sffamily \listofreviews}
		
		\cleardoublepage
	}
	
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Notation conventions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	\fancyhead[LE,RO]{\sffamily Notation}
	
	\input{notation}
	
	\cleardoublepage
	
	\pagenumbering{arabic}\setcounter{page}{1}	% For the remainder use arabic numeration
	
	\fancyhead[LE]{\sffamily \leftmark}					% Chapter in the left on even pages
	\fancyhead[RO]{\sffamily \nouppercase{\rightmark}}	% Section in the right on odd pages
	\pagestyle{fancy}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Bibliography
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newcommand{\makebibliography}[1]{
	\fancyhead[LE,RO]{\sffamily \ifenglish{Bibliography}{Literaturverzeichnis}}
	
	\phantomsection
	\ifenglish{
		\addcontentsline{toc}{chapter}{Bibliography}
	}{	
		\addcontentsline{toc}{chapter}{Literaturverzeichnis}
	}
	
	\ifenglish{
		\bibliographystyle{plain}
	}{
		\bibliographystyle{geralpha}
	}
	
	\bibliography{#1}
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Footnotes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newlength  {\footnotewidth}
\setlength  {\footnotewidth}{\textwidth}
\addtolength{\footnotewidth}{-5mm}

\makeatletter

\renewcommand{\@makefntext}[1]{
	\noindent
	\parbox[t]{5mm}{\mbox{\@thefnmark}}\parbox[t]{\footnotewidth}{#1}
}

\makeatother


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%           Misc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\makeatletter

\def\@makechapterhead#1{
  \vspace*{50\p@}
  {\parindent \z@ \raggedright \sffamily
        {\Large \MakeUppercase{\@chapapp} \space \thechapter}
        \par\nobreak
        \vskip 10\p@
    \interlinepenalty\@M
    \huge\bfseries #1\par\nobreak
    \vskip 40\p@
  }}

\def\@makeschapterhead#1{
  \vspace*{50\p@}
  {\parindent \z@ \raggedright
    \sffamily
    \interlinepenalty\@M
    \Huge \bfseries  #1\par\nobreak
    \vskip 40\p@
  }}

\renewcommand\section{\@startsection{section}{1}{\z@}
                                   {-3.5ex \@plus -1ex \@minus -.2ex}
                                   {2.3ex \@plus.2ex}
                                   {\sffamily\Large\bfseries}}
\renewcommand\subsection{\@startsection{subsection}{2}{\z@}
                                   {-3.25ex\@plus -1ex \@minus -.2ex}
                                   {1.5ex \@plus .2ex}
                                   {\sffamily\large\bfseries}}
\renewcommand\subsubsection{\@startsection{subsubsection}{3}{\z@}
                                   {-3.25ex\@plus -1ex \@minus -.2ex}
                                   {1.5ex \@plus .2ex}
                                   {\sffamily\normalsize\bfseries}}
\renewcommand\paragraph{\@startsection{paragraph}{4}{\z@}
                                   {3.25ex \@plus1ex \@minus.2ex}
                                   {-1em}
                                   {\sffamily\normalsize\bfseries}}
\renewcommand\subparagraph{\@startsection{subparagraph}{5}{\parindent}
                                   {3.25ex \@plus1ex \@minus .2ex}
                                   {-1em}
                                   {\sffamily\normalsize\bfseries}}
                                   
\makeatother

