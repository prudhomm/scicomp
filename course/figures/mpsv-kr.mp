%
% fig2dev (version 3.2.4) -L (m)mp version 0.05 --- Preamble
%

%
% mp output driver options:
%  -L mp -m 1.000000
%

% +MP-ADDITIONAL-HEADER
verbatimtex
%&latex
\documentclass{article}
\begin{document}
etex
% -MP-ADDITIONAL-HEADER


%SetFigFont macros for latex
verbatimtex
%
\begingroup\makeatletter\ifx\SetFigFont\undefined%
\gdef\SetFigFont#1#2#3#4#5{%
  \reset@font\fontsize{#1}{#2pt}%
  \fontfamily{#3}\fontseries{#4}\fontshape{#5}%
  \selectfont}%
\fi\endgroup%
\ifx\SetFigFontSize\undefined%
\gdef\SetFigFontSize#1#2{%
  \fontsize{#1}{#2pt}%
  \selectfont}%
\fi%
etex

% Make arrowheads mitered by default
% NOTE: subject to change (edited from plain.mp)
    def forwarr(text t) expr p =
      _apth:=p;_finarrf(t)
    enddef;
    def backarr(text t) expr p =
      _apth:=p;_finarrb(t)
    enddef;
    def _finarrf(text s) text t =
      if (s=0):fill arrowhead _apth  t withcolor white
      else: fill arrowhead _apth  t fi;
      linejoin:=0;
      draw arrowhead _apth  t
    enddef;
    def _finarrb(text s) text t =
      if (s=0):fill arrowhead reverse _apth  t withcolor white
      else: fill arrowhead reverse _apth  t fi;
      linejoin:=0;
      draw arrowhead reverse _apth  t
    enddef;


% Now draw the figure
beginfig(0)
% Some reasonable defaults
  ahlength:=7;
  ahangle:=30;
  labeloffset:=0;
  truecorners:=1;
  bboxmargin:=0;
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.56, 706.96)
    --(297.64, 706.96)
    --(297.64, 650.27)
    --(127.56, 650.27)
    --(127.56, 706.96)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.56, 650.27)
    --(297.64, 650.27)
    --(297.64, 593.57)
    --(127.56, 593.57)
    --(127.56, 650.27)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.56, 593.57)
    --(297.64, 593.57)
    --(297.64, 536.88)
    --(127.56, 536.88)
    --(127.56, 593.57)--cycle;
  draw p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (297.64, 536.88)
    --(368.50, 536.88);
  draw p withcolor (0.00,0.00,0.00);
  ahlength:=7.56;
  ahangle:=angle (7.56,1.89) * 2.0;
  forwarr(0) p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin polyline object
  linecap:=0;
  linejoin:=0;
  pickup pencircle scaled 0.94;
  path p;
  p = (127.56, 706.96)
    --(127.56, 763.65);
  draw p withcolor (0.00,0.00,0.00);
  ahlength:=7.56;
  ahangle:=angle (7.56,1.89) * 2.0;
  forwarr(0) p withcolor (0.00,0.00,0.00);
% End polyline object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{26}{31.3}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt $x$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((354.33,508.54))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{26}{31.3}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt $\hat{k}_3(y)$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((184.25,551.06))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{26}{31.3}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt $\hat{k}_2(y)$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((184.25,607.75))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{26}{31.3}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt $\hat{k}_1(y)$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((184.25,678.61))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
% Begin text object
  picture q;
  verbatimtex
   \def\mpsetfnt{%
       \SetFigFont{26}{31.3}{\familydefault}{\mddefault}{\updefault}%
    }%
  etex;
  q = btex \mpsetfnt $y$ etex;
  picture p;
  p = q rotated 0.00;
  label.urt(p,((85.04,749.48))+llcorner p) withcolor (0.00,0.00,0.00);
% End text object
endfig;
end
