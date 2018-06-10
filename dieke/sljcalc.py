from __future__ import division
import pdb

# Some code translated from Mike Reid's Pascal Code "linuxemp"
# Specifically zeeman operator stuff from sljcalc.p

from math import sqrt
from .wigner import Wigner6j


def odd(x):
    """
    returns False if x%2==0, otherwise returns true
    """
    if x % 2 == 0:
        return False
    else:
        return True


def istriad(j1, j2, j3):
    """
    Translated from Mike's emp code j1,j2,j3 integer
    Checks to see if j1 x j2 contains j3
    """
    if (j1+j2 >= j3) and (abs(j1-j2) <= j3) and (not odd(j1+j2+j3)):
        return True
    else:
        return False


def sixj(ja, jb, jc, jd, je, jf):
    """
    {computes 6j symbol arguments are twice actual value }
    { doesnt check validity }
    { rather than translating the pascall I just use the python version}
    """
    return Wigner6j(ja/2.0, jb/2.0, jc/2.0, jd/2.0, je/2.0, jf/2.0)


def t0kkval(k, rme, s1, l1, j1, s2, l2, j2):
    """
    {calculates slj reduced matrix element of uk  etc }
    """
    if (s1 != s2):
        return 0.0
    if not istriad(l1, l2, k):
        return 0.0
    if not istriad(j1, j2, k):
        return 0.0
    if odd((s2+l2+j1+k) // 2):
        sign = -1
    else:
        sign = 1
    return sign*sqrt((j1+1)*(j2+1))*sixj(j1,k,j2,l2,s1,l1)*rme;


def tk0kval(k, rme, s1, l1, j1, s2, l2, j2):
    """
    {calculates slj reduced matrix element of s   etc }
    """
    if not istriad(s1, s2, k):
        return 0.0
    if (l1 != l2):
        return 0.0
    if not istriad(j1, j2, k):
        return 0.0
    if odd((s1+l1+j2+k) // 2):
        sign = -1
    else:
        sign = 1
    return sign*sqrt((j1+1)*(j2+1))*sixj(j1,k,j2,s2,l1,s1)*rme
    

def tmagmomval(s1,l1,j1,s2,l2,j2):
    """
    calculates slj reduced matrix element of l+gs*s 
    uses formula for <]]s[[> and <]]l[[> and the 
    tk0k and t0kk functions
    """
    gs = 2.00232
    zero = False
    if (s1 != s2):
        return 0.0
    if (l1 != l2):
        return 0.0
    if (not istriad(j1, j2, 2)):
        return 0.0
    sx = s1/2.0
    lx = l1/2.0
    tmagmomval = gs * tk0kval(2,sqrt(sx*(sx+1)*(2*sx+1)),s1,l1,j1,s2,l2,j2)+t0kkval(2,sqrt(lx*(lx+1)*(2*lx+1)),s1,l1,j1,s2,l2,j2)


def reducedL(s1, l1, j1, s2, l2, j2):
    """
    Calculates the slj reduced matrix elements for L
    guessed from Mike's tmagmomval function.
    """
    if (s1 != s2):
        return 0.0
    if (l1 != l2):
        return 0.0
    if (not istriad(j1, j2, 2)):
        return 0.0
    # sx = s1/2.0
    lx = l1/2.0
    return t0kkval(2, sqrt(lx*(lx+1)*(2*lx+1)), s1, l1, j1, s2, l2, j2)


def reducedS(s1, l1, j1, s2, l2, j2):
    """
    Calculates the slj reduced matrix elements for S
    guessed from Mike's tmagmomval function.
    """
    if (s1 != s2):
        return 0.0
    if (l1 != l2):
        return 0.0
    if (not istriad(j1, j2, 2)):
        return 0.0
    sx = s1/2.0
    # lx = l1/2.0
    return tk0kval(2, sqrt(sx*(sx+1)*(2*sx+1)), s1, l1, j1, s2, l2, j2)


oldcode =  """
program sljcalc(input, {altinput,} output{, textin, textin2,
		  textout, diskmein, diskmeout});

{ uses readcross o/p to calculate slj (reduced) matrix elements }
{ deals with upper triangle of matrix only }
{ author: mike reid }

{ Nov 27 92 changes so3labeltype to integer in repsonse to an apparent
  bug in SUN Pascal }

{ Nov 14 92 - Major revision to readcr and sljcalc:
  the factor:
        (-1)^k sqrt(2*k+1)
  is removed from the calculation routine tkk0calc
      so that , e.g. zeta is (sl)0, not s.l.  }

{ Dec 6 94 (GWB) - Minor format modifications for compatability with
  new READCR two-configuration-allowed output, and variable-length
  state labels. }

{ November 17, 1999 (MFR) added procedures to read crosswhite SLJ files
 and write RACAH files. }
{ July 20 2000 (MFR) minor fixes to allow for large numbers of states }
{ Oct 12 2000 (MFR) added code to calculate trace of operators }
{ Dec 07 2000 (MFR) added ALL option to makecfitinput to allow free-ion fd calculations}
{ Oct 07 2001 (MFR) added special case so that W101 is calcualated
 using tsljval. Note that we can calculate S from sqrt(7/2)*W101.}

{DOLLARi tp.def}
{$ifndef turbo}
{$include "gnu.def"}
{$endif}

const
{DOLLARi nj.cns}
{$ifndef turbo}
{$include "gnu_io.cns"}
{$include "runner.cns"}
{$include "nj.cns"}
{$endif}
    progname = 'SLJCALC 07 Oct 01   '; {must be 20 char}

{turbo}
(*
    maxnsl  =  64;
    maxnslj = 128;
*)
{large model}
(*
    maxnsl  =  600;
    maxnslj = 1000;
*)
{huge model}
(**)
    maxnsl  = 1237;
    maxnslj = 4000;
(**)  
{lpretty label length}
    maxsllab  = 10;
    maxsljlab = 12;

type
{$ifndef turbo}
{$include "gnu_io.typ"}
{$include "runner.typ"}
{$endif}
    (* so3labeltype = 0..63; *)
    so3labeltype = integer; 
    alfasl  = packed array [1..maxsllab]  of char;
    alfaslj = packed array [1..maxsljlab] of char;
    sllabtype = (*packed*) record
                    sly, sls, sll, slx: so3labeltype;
                    slpretty: alfasl;
                  end;
    sljlabtype = (*packed*) record
                    sljy, sljs, sljl, sljx, sljj: so3labeltype;
                    sljwhichsl:1..maxnsl;
                    sljpretty: alfaslj;
                  end;
    jlabeltype = array[1..10] of integer;

var
{DOLLARi nj.var}
{$ifndef turbo}
{$include "gnu_io.var"}
{$include "runner.var"}
{$include "nj.var"}
{$endif}
    textin,textin2,textout: text;
    diskmein,diskmeout:rmefile;
    nelectrons,nsl, nslj,iuu,ivv,juu,jvv,item,maxrange:integer;
    sllab: array[1..maxnsl] of sllabtype;
    sljlab: array[1..maxnslj] of sljlabtype;
    redmatrix: array[1..maxnsl, 1..maxnsl] of real;
    sljmatrix: array[1..maxnslj, 1..maxnslj] of real;
    maden,twoconfig :boolean;
    lengthsl,lengthslj:integer;
    termlab: packed array [1..20] of char;

{ PROCEDURES }

{DOLLARi nj.zp}
{$ifndef turbo}
{$include "gnu_io.zp"}
{$include "runner.zp"}
{$include "nj.zp"}
{$endif}

procedure diagnost;
begin
  writeln;
end;

function tkk0val(k:integer; rme:real; s1,l1,j1,s2,l2,j2:integer):real;
{calculates matrix element [not reduced] of v11 etc }
{ Nov 14 92 divide by (-1)^k * sqrt(2*k+1) so that we are NOT 
  doing the calculation for a.b, but for (ak bk)0 }
var
  sign:integer;
  zero:boolean;
begin
  zero:=false;
  if j1<>j2 then zero:=true else
  if not istriad(s1,s2,k) then zero:=true else
  if not istriad(l1,l2,k) then zero:=true;
  if zero then tkk0val:=0 else begin
    if odd((j1+l1+s2-k) div 2) then sign:=-1 else sign:=1;
    tkk0val:=sign*sixj(l1,l2,k,s2,s1,j1) * rme / sqrt(k+1);
  end;
end {tkk0val};


function t0kkval(k:integer; rme:real; s1,l1,j1,s2,l2,j2:integer):real;
{calculates slj reduced matrix element of uk  etc }
var
  sign:integer;
  zero:boolean;
begin
  zero:=false;
  if s1<>s2 then zero:=true else
  if not istriad(l1,l2,k) then zero:=true else
  if not istriad(j1,j2,k) then zero:=true;
  if zero then t0kkval:=0 else begin
    if odd((s2+l2+j1+k) div 2) then sign:=-1 else sign:=1;
    t0kkval:=sign*sqrt((j1+1)*(j2+1))*sixj(j1,k,j2,l2,s1,l1)*rme;
  end;
end {t0kkval};

function tk0kval(k:integer; rme:real; s1,l1,j1,s2,l2,j2:integer):real;
{calculates slj reduced matrix element of s   etc }
var
  sign:integer;
  zero:boolean;
begin
  zero:=false;
  if not istriad(s1,s2,k) then zero:=true else
  if l1<>l2 then zero:=true else
  if not istriad(j1,j2,k) then zero:=true;
  if zero then tk0kval:=0 else begin
    if odd((s1+l1+j2+k) div 2) then sign:=-1 else sign:=1;
    tk0kval:=sign*sqrt((j1+1)*(j2+1))*sixj(j1,k,j2,s2,l1,s1)*rme;
  end;
end {tk0kval};

function tsljval(s,l,j:integer; rme:real; s1,l1,j1,s2,l2,j2:integer):real;
{calculates slj reduced matrix element of general tslj }
var
  zero:boolean;
begin
  zero:=false;
  if not istriad(s1,s2,s) then zero:=true else
  if not istriad(l1,l2,l) then zero:=true else
  if not istriad(j1,j2,j) then zero:=true;
  if zero then tsljval:=0 else begin
    tsljval:=sqrt((j1+1)*(j2+1)*(j+1))*rme
            *ninej(s1,s2,s,l1,l2,l,j1,j2,j);
  end;
end {tsljval};

function tmagmomval(s1,l1,j1,s2,l2,j2:integer):real;
{calculates slj reduced matrix element of l+gs*s }
{ use formula for <]]s[[> and <]]l[[> and the    }
{ above tk0k and t0kk functions }
const
  gs=2.00232;
var
  zero:boolean;
  sx,lx:real;
begin
  zero:=false;
  if s1<>s2 then zero:=true else
  if l1<>l2 then zero:=true else
  if not istriad(j1,j2,2) then zero:=true;
  if zero then tmagmomval:=0 else begin
    sx:=s1/2; lx:=l1/2;
    tmagmomval:= gs *
                 tk0kval(2,sqrt(sx*(sx+1)*(2*sx+1)),
                         s1,l1,j1,s2,l2,j2)    +
                 t0kkval(2,sqrt(lx*(lx+1)*(2*lx+1)),
                         s1,l1,j1,s2,l2,j2);
  end;
end {tmagmomval};

procedure openipandop;
{ open textin and textout and diskmeout for read/write }
var
  title, title2 :titletype;
  ar, i: integer;
begin
  if not gettitle(title) then warn('TITLE REQUIRED      ', 'openipando');
  tabout; write('READING FROM ');  writetitle(output, title); writeln;
  checkpresent(title);
  textreset(textin, title);
  if not gettitle(title) then warn('O/P TITLE REQUIRED  ', 'openipando');
  title2:=title;
  addsuffix(title, sljinfosuffix);
  tabout; write('WRITING INFO TO ');  writetitle(output, title); writeln;
  textrewrite(textout, title);
  addsuffix(title2, mebinsuffix);
  tabout; write('WRITING MES  TO ');  writetitle(output, title2); writeln;
  rmerewrite(diskmeout, title2);
end {openipandop};

procedure findmaxrange; { return then max of (iuu,ivv,juu,jvv) }
var
  mr:integer;
begin
  mr:=0;
  if iuu>mr then mr:=iuu;
  if ivv>mr then mr:=ivv;
  if juu>mr then mr:=juu;
  if jvv>mr then mr:=jvv;
  maxrange:=mr;
end {maxrange};

procedure readsllabs;
{ read sl labels, write to textout and output }
var
  wlab,bcol,i:integer;
  buff:cardimage;
begin
  twoconfig:=false;
  tabout; write('READING STATES'); writeln;
  getfilecard(textin,buff); tabout; write(buff);writeln; bcol:=1;
  nelectrons:=readint(buff,bcol);
  getfilecard(textin,buff); tabout; write(buff);writeln;
  getfilecard(textin,buff); bcol:=1;
  nsl:=readint(buff,bcol); 
  if nsl>maxnsl then warn('increase maxnsl     ', 'readsllabs');
  skipblanks(buff,bcol); skiptoblank(buff,bcol);
  lengthsl:=readint(buff,bcol);
  if lengthsl>maxsllab then warn('increase maxsllab   ','readsllabs');
  writeln; tabout;
  write(nsl:5,' STATES');
  for i:=7 to maxsllab do write(' ');
  write('   Y  2S  2L   X'); writeln; writeln;
  maden:=true;
  for wlab:=1 to nsl do with sllab[wlab] do begin
    getfilecard(textin,buff);
    for i:=1 to maxsllab do slpretty[i]:=buff[i];
    bcol:=maxsllab+1;
    sly:=readint(buff,bcol);
     if sly<>0 then twoconfig:=true;
    sls:=readint(buff,bcol);
    sll:=readint(buff,bcol);
    slx:=readint(buff,bcol);
    if (maxrange=0)or(wlab<=maxrange) then begin
      tabout; write(wlab:5,' ',slpretty,sly:4,sls:4,sll:4,slx:4);writeln
    end else
      if maden then begin
        for i:=1 to 3 do begin tabout; write('      .'); writeln; end;
        maden:=false;
      end;
  end;
  prstatus(output);
end {readsllabs};

procedure sljlabgen;
{ generate slj labels. write to textout and output }
var
  wsl,wslj,j,i,jlab,jmax:integer;
begin
  writeln; tabout; write('SLJ LABELS'); writeln; writeln;
  wslj:=0;
  for wsl:=1 to nsl do
    with sllab[wsl] do begin
      j:=abs(sls-sll); jmax:=sls+sll;
      while j<=jmax do begin
        wslj:=wslj+1;
        if wslj>maxnslj then warn('increase maxnslj    ', 'sllabgen  ');
        with sljlab[wslj] do begin
          sljy:=sly; sljs:=sls; sljl:=sll; sljx:=slx; sljj:=j;
          sljwhichsl:=wsl;
          for i:=1 to maxsllab do sljpretty[i]:=slpretty[i];
          for i:=maxsllab+1 to maxsljlab do sljpretty[i]:=' ';
          if odd(j) then jlab:=j else jlab:=j div 2;
          lengthslj:=lengthsl+2;
          sljpretty[lengthslj]:=sljpretty[lengthsl];
          sljpretty[lengthsl]:=' ';    
          if jlab>=10 then
            sljpretty[lengthslj-2]:=chr((jlab div 10)+ord('0'));
          sljpretty[lengthslj-1]:=chr((jlab mod 10)+ord('0'));
        end;
        j:=j+2;
      end;
    end;
  write(textout,nelectrons:3,' ELECTRONS: SLJ STATES, MES AND RMES');
  writeln(textout);
  writecreationdate(textout, progname);
  nslj:=wslj;
  tabout; write(nslj:5,' STATES');
  for i:=7 to maxsljlab do write(' ');
  write('   Y  2S  2L   X  2J WHICHSL'); writeln; writeln;
  write(textout,nslj:4,'  STATES, ',lengthslj:4,'  CHARACTER LABELS'); 
  writeln(textout);
  maden:=true;
  for wslj:=1 to nslj do
    with sljlab[wslj] do begin
      if (maxrange=0)or(wslj<=maxrange) then begin
        tabout; write(wslj:5,' ',sljpretty,sljy:4,sljs:4,sljl:4,sljx:4,sljj:4
          ,sljwhichsl:4); writeln;
      end else
        if maden then begin
          for i:=1 to 3 do begin tabout; writeln('      .');end;
          maden:=false;
        end;
      write(textout,sljpretty, sljy:6,sljs:6,sljl:6,sljx:6,sljj:6
        ,sljwhichsl:6); writeln(textout)
    end;
  writeln;
  prstatus(output);
end {sljlabgen};

procedure readrme(var tenlab:alfa; var ts,tl,ntj:integer;
                  var tj:jlabeltype);
{ read reduced matrix elements of a tensor from textin }
var
  buff,buff1:cardimage;
  bcol,wme,nslme,i,j:integer;
  r:real;
begin
  for i:=1 to nsl do for j:=i to nsl do redmatrix[i,j]:=0;
  getfilecard(textin,buff);
  bcol:=1; getstr(buff,bcol,buff1);
  for i:=1 to 10 do tenlab[i]:=buff1[i];
  ts:=readint(buff,bcol);
  tl:=readint(buff,bcol);
  skipblanks(buff,bcol); ntj:=0;
  while bcol<=cardwidth do begin
    ntj:=ntj+1;
    tj[ntj]:=readint(buff,bcol); skipblanks(buff,bcol);
  end;
  writeln; writeln; writeln;
  tabout; write('*****  ');
  write(tenlab, 'S=', ts div 2:2, '    L=',tl div 2:2, '   J=');
  for i:=1 to ntj do begin
    write(tj[i] div 2 :2); if i<>ntj then write(',');
  end;
  write('    *****'); writeln;
  getfilecard(textin,buff); bcol:=1;
  nslme:=readint(buff,bcol);
  tabout; write(nslme:5,' SL (REDUCED) MATRIX ELEMENTS'); writeln;
  for wme:=1 to nslme do begin
    read(textin,i,j,r);
    redmatrix[i,j]:=r;
    (*if debug then writeln(i:5, j:5, r:12);*)
  end;
  if nslme<>0 then readln(textin);
end {readrme};

function inrange(i,j:integer):boolean;
begin
  inrange:=((iuu<=i)and(i<=ivv)and(juu<=j)and(j<=jvv))or(iuu*ivv*juu*jvv=0)
end {inrange};

procedure matelcalc;
{ calculate matrix elements }
var
  sljlabel				       : alfaslj;
  tenlabel,lab				       : alfa;
  ts,tl,ntj,wtj,wsljme,wslja,wsljb,wsla,wslb,i : integer;
  tcol,jlabel				       : integer;
  tj					       : jlabeltype;
  mel,rme,ss				       : real;
  sljme					       : rmetype;
  trace					       : real; 


procedure debugwrite(sljme: rmetype; i:integer);
      var
        k:integer;
      begin
        with sljme do
          if inrange(mi,mj) then begin
            item:=item+1; if item=1 then tabout;
            sljlabel:=sljlab[mi].sljpretty; 
            sljlabel[1]:='<'; sljlabel[lengthslj]:=']';
            for k:=1 to lengthslj do write(sljlabel[k]);
            for k:=1 to lengthslj do write(sljlab[mj].sljpretty[k]);
            write('=',me:12:6,'   ');
          { write(sljlabel,sljlab[mj].sljpretty,'=',me:12:6,'   '); }
            if (item mod 3)=0 then begin writeln; tabout; end;
          end;
      end { debugwrite };


begin { matelcalc }
  tabout; write('MATRIX ELEMENTS'); writeln; writeln;
  tabout;
  write('*****  EAVG     ', nslj:5, ' MATRIX ELEMENTS *****');
  writeln;writeln;
  write(textout, 'EAVG         0   0   0'); writeln(textout);
  write(textout, nslj:5);  writeln(textout);
  item:=0;
  for i:=1 to nslj do begin
    with sljme do begin
      mi:=i; mj:=i; me:=1;
    end;
    rmewrite(diskmeout, sljme);
    (*if debug then debugwrite(sljme, i);*)
  end;
  if (debug)and((item mod 3)<>0) then diagnost;
  if twoconfig then begin
    item:=0;
    for i:=1 to nslj do
     if sljlab[i].sljy<>0 then item:=item+1;
    tabout;
    write('*****  EEXC     ', item:5, ' MATRIX ELEMENTS *****');
    writeln;
    write(textout, 'EEXC         0   0   0'); writeln(textout);
    write(textout, item:5);  writeln(textout);
    item:=0;
    for i:=1 to nslj do 
    if sljlab[i].sljy<>0 then begin
      with sljme do begin
        mi:=i; mj:=i; me:=1;
      end;
      rmewrite(diskmeout, sljme);
      if debug then debugwrite(sljme, i);
    end;
    if (debug)and((item mod 3)<>0) then diagnost;
  end;
  while not eof(textin) do begin
    readrme(tenlabel,ts,tl,ntj,tj);
    trace:=0;
    if ntj<>1 then begin
      tcol:=1;
      while (tenlabel[tcol]<>' ') and (tcol<9) do tcol:=tcol+1;
    end;
    for wtj:=1 to ntj do begin
      if ntj<>1 then begin
        for i:=tcol to 10 do tenlabel[i]:=' ';
        jlabel:=tj[wtj] div 2; i:=tcol;
        if jlabel>10 then begin
          tenlabel[i]:=chr(jlabel div 10 + ord('0') );
          i:=i+1;
        end;
        tenlabel[i]:=chr(jlabel mod 10 + ord('0') );
      end;
      write(textout,tenlabel,ts:4,tl:4,tj[wtj]:4); writeln(textout);

      lab:=tenlabel; for i:=5 to 10 do lab[i]:=' ';

{ sccf  put in delta(s1,s2) * sqrt(s(s+1)/(2s+1)) factor }
{       nb: transforms as s=0, l=k, j=k                 }

      if lab = 'SCCF      ' then begin
        tabout;
        write('MULTIPLYING V1K RME BY DELTA(S1,S2)*SQRT(S(S+1)/(2S+1))');
        writeln;
        wsljme:=0; item:=0;
        for wslja:=1 to nslj do
          for wsljb:=wslja to nslj do
            if (sljlab[wslja].sljs=sljlab[wsljb].sljs) {delta(s1,s2)} then
              if (sljlab[wslja].sljs<>0) then begin
                wsla:=sljlab[wslja].sljwhichsl;
                wslb:=sljlab[wsljb].sljwhichsl;
                rme:=redmatrix[wsla,wslb];
                if rme<>0 then begin
                  ss:=sljlab[wslja].sljs /2;
                  rme:=rme * sqrt(ss*(ss+1)/(2*ss+1)); { v1k->sccf factor }
                  mel:=t0kkval(tl,rme,
                    sljlab[wslja].sljs,sljlab[wslja].sljl,sljlab[wslja].sljj,
                    sljlab[wsljb].sljs,sljlab[wsljb].sljl,sljlab[wsljb].sljj);
                  if mel<>0 then begin
                    wsljme:=wsljme+1;
                    with sljme do begin
                      mi:=wslja; mj:=wsljb; me:=mel;
		      if mi=mj then trace:=trace+mel* sqrt(sljlab[wslja].sljj+1);
                    end;
                    rmewrite(diskmeout, sljme);
                    if debug then debugwrite(sljme, wsljme);
                  end {if mel};
                end {if rme};
              end {if (sljlab...};
        if (debug)and((item mod 3)<>0) then diagnost;
      end else

{ magmoment }

      if tenlabel = 'MAGMOM    ' then begin
        tabout; write('USING FORMULAE FOR MAGNETIC MOMENT'); writeln;
        wsljme:=0; item:=0;
        for wslja:=1 to nslj do
         for wsljb:=wslja to nslj do
          if sljlab[wslja].sljwhichsl = sljlab[wsljb].sljwhichsl then
           { must have same parentage ! }
           begin
             wsla:=sljlab[wslja].sljwhichsl;
             wslb:=sljlab[wsljb].sljwhichsl;
             mel:=tmagmomval(
                sljlab[wslja].sljs,sljlab[wslja].sljl,sljlab[wslja].sljj,
                sljlab[wsljb].sljs,sljlab[wsljb].sljl,sljlab[wsljb].sljj);
             if mel<>0 then begin
               wsljme:=wsljme+1;
               with sljme do begin
                 mi:=wslja; mj:=wsljb; me:=mel;
		 if mi=mj then trace:=trace+mel* sqrt(sljlab[wslja].sljj+1);
               end;
              rmewrite(diskmeout, sljme);
              if debug then debugwrite(sljme, wsljme);
             end;
           end;
        if (debug)and((item mod 3)<>0) then diagnost;
      end else

{ sl scalars }

      if (ts=0) and (tl=0) and (tj[wtj]=0) then begin
	 if debug then writeln('sl scalar');
        wsljme:=0; item:=0;
        for wslja:=1 to nslj do
        for wsljb:=wslja to nslj do
        begin
          if sljlab[wslja].sljj = sljlab[wsljb].sljj then begin
            wsla:=sljlab[wslja].sljwhichsl;
            wslb:=sljlab[wsljb].sljwhichsl;
            rme:=redmatrix[wsla,wslb];
            if rme<>0 then begin
              wsljme:=wsljme+1;
              with sljme do begin
                mi:=wslja; mj:=wsljb; me:=rme;
		if mi=mj then begin
		  trace:=trace+me * (sljlab[wslja].sljj+1);
		  (*
		  if debug then begin
		    writeln; 
		    writeln(mi:6, me:12:4, ' trace:', trace:12:4);
		  end;
		  *) 
		end;
              end;
              rmewrite(diskmeout, sljme);
              if debug then debugwrite(sljme, wsljme);
            end;
          end;
        end;
        if (debug)and((item mod 3)<>0) then diagnost;
      end else

{ tkk0 }

      if tj[wtj]=0 then begin
        wsljme:=0; item:=0;
        for wslja:=1 to nslj do
        for wsljb:=wslja to nslj do
        begin
          wsla:=sljlab[wslja].sljwhichsl;
          wslb:=sljlab[wsljb].sljwhichsl;
          rme:=redmatrix[wsla,wslb];
          if rme<>0 then begin
            mel:=tkk0val(ts,rme,
               sljlab[wslja].sljs,sljlab[wslja].sljl,sljlab[wslja].sljj,
               sljlab[wsljb].sljs,sljlab[wsljb].sljl,sljlab[wsljb].sljj);
            if mel<>0 then begin
              wsljme:=wsljme+1;
              with sljme do begin
                mi:=wslja; mj:=wsljb; me:=mel;
		if mi=mj then trace:=trace+mel * (sljlab[wslja].sljj+1);
              end;
              rmewrite(diskmeout, sljme);
              if debug then debugwrite(sljme, wsljme);
            end;
          end;
        end;
        if (debug)and((item mod 3)<>0) then diagnost;
      end else

{ t0kk }

      if ts=0 then begin
        wsljme:=0; item:=0;
        for wslja:=1 to nslj do
        for wsljb:=wslja to nslj do
        begin
          wsla:=sljlab[wslja].sljwhichsl;
          wslb:=sljlab[wsljb].sljwhichsl;
          rme:=redmatrix[wsla,wslb];
          if rme<>0 then begin
            mel:=t0kkval(tl,rme,
              sljlab[wslja].sljs,sljlab[wslja].sljl,sljlab[wslja].sljj,
              sljlab[wsljb].sljs,sljlab[wsljb].sljl,sljlab[wsljb].sljj);
            if mel<>0 then begin
              wsljme:=wsljme+1;
              with sljme do begin
                mi:=wslja; mj:=wsljb; me:=mel;
		if mi=mj then trace:=trace+mel* sqrt(sljlab[wslja].sljj+1);
              end;
              rmewrite(diskmeout, sljme);
              if debug then debugwrite(sljme, wsljme);
            end;
          end;
        end;
        if (debug)and((item mod 3)<>0) then diagnost;
      end else

{ tk0k } {note that W101 is explicity excluded here }

      if (tl=0) and (not (tenlabel = 'W101      ')) then begin
        wsljme:=0; item:=0;
        for wslja:=1 to nslj do
        for wsljb:=wslja to nslj do
        begin
          wsla:=sljlab[wslja].sljwhichsl;
          wslb:=sljlab[wsljb].sljwhichsl;
          rme:=redmatrix[wsla,wslb];
          if rme<>0 then begin
            mel:=tk0kval(ts,rme,
              sljlab[wslja].sljs,sljlab[wslja].sljl,sljlab[wslja].sljj,
              sljlab[wsljb].sljs,sljlab[wsljb].sljl,sljlab[wsljb].sljj);
            if mel<>0 then begin
              wsljme:=wsljme+1;
              with sljme do begin
                mi:=wslja; mj:=wsljb; me:=mel;
		if mi=mj then trace:=trace+mel* sqrt(sljlab[wslja].sljj+1);
              end;
              rmewrite(diskmeout, sljme);
              if debug then debugwrite(sljme, wsljme);
            end;
          end;
        end;
        if (debug)and((item mod 3)<>0) then diagnost;
      end else

{ tslj }

      begin
	if tenlabel = 'W101      ' then begin
          tabout; write('USING TSLJVAL TO CALCULATE W101'); writeln;
        end;  
        wsljme:=0; item:=0;
        for wslja:=1 to nslj do
        for wsljb:=wslja to nslj do
        begin
          wsla:=sljlab[wslja].sljwhichsl;
          wslb:=sljlab[wsljb].sljwhichsl;
          rme:=redmatrix[wsla,wslb];
          if rme<>0 then begin
            mel:=tsljval(ts,tl,tj[wtj],rme,
               sljlab[wslja].sljs,sljlab[wslja].sljl,sljlab[wslja].sljj,
               sljlab[wsljb].sljs,sljlab[wsljb].sljl,sljlab[wsljb].sljj);
            if mel<>0 then begin
              wsljme:=wsljme+1;
              with sljme do begin
                mi:=wslja; mj:=wsljb; me:=mel;
		if mi=mj then trace:=trace+mel* sqrt(sljlab[wslja].sljj+1);
              end;
              rmewrite(diskmeout, sljme);
              if debug then debugwrite(sljme, wsljme);
            end;
          end;
        end;
        if (debug)and((item mod 3)<>0) then diagnost;
      end;

      tabout;
      write('*****  ',tenlabel,wsljme:5,' MATRIX ELEMENTS  *****');
      write(' TRACE: ',trace:6:3,' *****');
      writeln;
      write(textout,wsljme:5); writeln(textout);
      prstatus(output);

    end { for wtj };
  end { while not eof };
  rmeclose(diskmeout);
  textclose(textin);
  textclose(textout);
end {matelcalc};

procedure sljcalculation;
{ orchestrate the calculation }
begin
  openipandop;
  readsllabs;
  sljlabgen;
  matelcalc;
end {sljcalculation};

procedure makecfitinput;
{ prepare output for cfit }
{ read from o/p of sljcalculation write .st, .m$i }
{ use same me file as before, ie .m$m already created }
{ assumes all free-ion tensors are before all others and   }
{ stops when a non-free-ion tensor is encountered          }
var
  titlei,titleout, titleout2:titletype;
  bcol,i,j,nme,ts,tl,tj,ar:integer;
  buff,buff1:cardimage;
  lab:alfa;
  alltensors,finished:boolean;
begin
  if not gettitle(titlei) then warn('TITLE REQUIRED      ', 'makecfitin');
  titleout:=titlei; titleout2:=titlei;
  addsuffix(titlei, sljinfosuffix);
  tabout; write('READING INFO FROM ');  writetitle(output, titlei);writeln;
  addsuffix(titleout, statesuffix);
  tabout; write('WRITING STATES TO ');  writetitle(output, titleout);writeln;
  textrewrite(textout, titleout);
  checkpresent(titlei);
  textreset(textin, titlei); { done in this order to avoid file in use problems }

  alltensors:=false;
  nextword;
  while incol<cardwidth do begin
    if samestring(inbuff,incol,'ALL       ',4) then alltensors:=true else
    warn('INVALID OPTION      ', 'MAKECFITIN');
    nextword;
  end;
  
{ write states }

  getfilecard(textin,buff); tabout; write(buff); writeln;
  bcol:=1; nelectrons:=readint(buff,bcol);
  write(textout,nelectrons:3,' ELECTRONS. SLJ STATES'); writeln(textout);
  getfilecard(textin,buff); tabout; write(buff); writeln;
  writecreationdate(textout, progname);
  getfilecard(textin,buff); tabout; write(buff); writeln;
  write(textout,buff); writeln(textout);
  bcol:=1; nsl:=readint(buff,bcol);
  for i:=1 to nsl do begin;
    getfilecard(textin,buff);
    j:=cardwidth; while buff[j]=' ' do j:=j-1;
    write(textout,buff:j+1); writeln(textout);
  end;
  textclose(textout);

{ write matrix elements }
  addsuffix(titleout2, meinfosuffix);
  tabout; write('WRITING ME INFO TO   ');  writetitle(output, titleout2);
  writeln;
  textrewrite(textout, titleout2);
  write(textout,nelectrons:3,' ELECTRONS. SLJ MATRIX ELEMENTS');
  writeln(textout);
  writecreationdate(textout, progname);
  finished :=false;
  while not finished do begin
    getfilecard(textin,buff); bcol:=1;
    getstr(buff,bcol,buff1);
    for i:=1 to 10 do lab[i]:=buff1[i];
    ts:=readint(buff,bcol);
    tl:=readint(buff,bcol);
    tj:=readint(buff,bcol);
    if ((tj<>0) and (not alltensors)) then begin
      finished:=true; writeln; tabout;
      write('END OF FREE-ION TENSORS: FINISHING.'); writeln;
    end else begin
      tabout; write(buff); writeln;
      getfilecard(textin,buff); bcol:=1; nme:=readint(buff,bcol);
      write(textout,lab,' ',nme:5); writeln(textout);
    end;
    if eof(textin) then finished:=true;
  end;
  textclose(textin);
  textclose(textout);
  prstatus(output);
end {makecfitinput};

procedure setrange;
begin
  iuu:=nextint; ivv:=nextint; juu:=nextint; jvv:=nextint;
  findmaxrange;
end;


{ Procedures for replacing SLJ free-ion matrix elements by
  newer Crosswhite ones }

procedure setuptermlabs;
begin
  termlab:='SPDFGHIKLMNOQRT     ';
end {setuptermlabs};

procedure translatelab(var lab:alfa; var s, l, x:integer);
{ tranlates lab into s,l,x - s&l 2x actual value }
var i:integer;
begin
    s:=ord(lab[1])-ord('0')-1;
    if lab[3]=' ' then x:=0 else begin
      x:=ord(lab[3])-ord('0');
      if x=0 then x:=10;
    end;
    i:=1;
    while (lab[2]<>termlab[i]) and (termlab[i]<>' ') do i:=i+1;
    if termlab[i]=' ' then begin
      tabout; write(lab); warn('INVALID LABEL       ', 'translatel');
    end else l:=(i-1) * 2;
end{translatelab};

function labpos(s, l, x, j:integer): integer;
{ find a label in sljlab and return its position }
var
  low,high,place,y:integer;
begin
  y:=0; {only search in the fn section of the labels}
  low:=1; high:=nslj; place:=(low+high) div 2;
  repeat
    with sljlab[place] do begin
      if sljy>y then high:=place-1 else
      if sljy<y then low:=place+1 else
      if sljs>s then low:=place +1 else
      if sljs<s then high:=place-1 else
      if sljl>l then high:=place-1 else
      if sljl<l then low:=place +1 else
      if sljx>x then high:=place-1 else
      if sljx<x then low:=place +1 else
      if sljj>j then high:=place-1 else
      if sljj<j then low:=place +1 else
        begin high:=place; low:=place; end;
    end;
    place:=(low+high) div 2;
  until low>=high;
  with sljlab[place] do begin
    if (y=sljy) and (s=sljs) and (l=sljl) and (x=sljx) 
       and (j=sljj) then labpos:=place else
      begin
        tabout; write(y:3,s:3,l:3,x:3,j:3); 
        warn('CANNOT FIND LABEL   ', 'labpos    ');
      end;
  end;
end {labpos};

function getmatrix(var searchlab:alfa; var textin:text):boolean;
{ read the new crosswhite file and create a matrix in sljmatrix }
{ note that the f12 file had to be hand-edited to add the proper labels!}
var
  buff,buff1:cardimage; 
  lab,fflab:alfa;
  s, l, x, j, i, bcol, bcol1, njstates, wjstate, wsljstate: integer; 
  ji, jj, statei, statej: integer; 
  value: real; 
  sljtranslation: array[1..maxnslj] of integer; 
begin
  getmatrix:=false; 
  for statei := 1 to nslj do 
    for statej := 1 to nslj do
      sljmatrix[statei,statej]:=0;
  while not eof(textin) do begin
    getfilecard(textin,buff); bcol:=1; 
    (* tabout; write(buff); writeln;  *)
    { header }
    extractstr(buff,bcol,6,buff1); bcol1:=1; j:=readint(buff1,bcol1);
    if odd(nelectrons) then { j should be twice actual value }
       j := (j-1)*2 + 1
    else
       j := (j-1)*2;
    extractstr(buff,bcol,6,buff1); bcol1:=1; njstates:=readint(buff1,bcol1);
    (* writeln(j:3, njstates:3); *)
     for wjstate := 1 to njstates do begin
       getfilecard(textin,buff); bcol:=1; 
       (* tabout; write(buff); writeln; *)
       extractstr(buff,bcol,6,buff1); 
       extractstr(buff,bcol,6,buff1); 
       extractstr(buff,bcol,6,buff1); 
       extractstr(buff,bcol,6,buff1); 
       extractstr(buff,bcol,12,buff1); 
       extractstr(buff,bcol,2,buff1);
       extractstr(buff,bcol,3,buff1); 
       for i:=1 to 10 do lab[i]:=buff1[i];
       translatelab(lab, s, l, x); 
       wsljstate:=labpos(s, l, x, j); 
       sljtranslation[wjstate]:=wsljstate;
       (* writeln('state: ', wjstate:3, wsljstate:3, 
                   ' 2j= ', j:1, '  ', lab); *)
     end;
     repeat
       getfilecard(textin,buff); bcol:=1;  
       if not samestring(buff,bcol,'          ',10) then begin
         (* tabout; write(buff); writeln; *)
         extractstr(buff,bcol,6,buff1); 
         extractstr(buff,bcol,6,buff1); bcol1:=1; ji:=readint(buff1,bcol1);
         extractstr(buff,bcol,6,buff1); bcol1:=1; jj:=readint(buff1,bcol1);
         extractstr(buff,bcol,6,buff1); 
         extractstr(buff,bcol,12,buff1); bcol1:=1; value:=readnum(buff1,bcol1);
         extractstr(buff,bcol,1,buff1); 
         extractstr(buff,bcol,10,buff1); 
         for i:=1 to 10 do lab[i]:=buff1[i];
         if lab='.01ALPH   ' then begin
           lab:='ALPHA     ';
           value:=value*1000;
         end; 
         if lab='  T2    T3' then {some n>7 labels are wrong}
           lab:='T2        ';
         fflab:=lab;
         i:=1; 
         while (i<8) and (fflab[i]<>' ') do i:=i+1;
         fflab[i]:='F'; fflab[i+1]:='F';
         (* writeln(ji:3, jj:3, sljtranslation[ji]:3, sljtranslation[jj]:3,
                  value:15:8, ' ',  lab, fflab);  *)
         if (searchlab=lab) or (searchlab=fflab) then begin
           sljmatrix[sljtranslation[ji], sljtranslation[jj]]:=value; 
           getmatrix := true; 
         end; 
       end;
     until (samestring(buff,1,'          ',10) or eof(textin));	
  end;
end {getmatrix};


procedure ncrread;
{ read from sljcalc output and new crosswhite input and substitute the
  operators found in the new file. 
  syntax: ncrread <new N> <old slj file> <new crosswhite file> <new slj file>
  note that the slj "files" are actually three files 
  Also note that the signs of Uk and V1k operators are NOT changed!
  Only important change is T2 for N=6 to 12 (yes, 12).
 Crosswhite files are called fNmp.dat. Had to hand-edit f12mp.dat to make it
 work. 
}
  var
    title, title2, crosstitle: titletype; 
    wslj, i, bcol: integer;
    buff, buff1: cardimage;
    finished : boolean;
    ts, tl, tj, nme, newnme, wme: integer;    
    lab : alfa; 
    sljme: rmetype;
    statei, statej: integer;  
    replacematrix: boolean;
    newnelectrons: integer; 

begin
  setuptermlabs;

  newnelectrons := nextint; 
  tabout; writeln('CREATING FILES FOR N = ', newnelectrons:1);

  if not gettitle(title) then warn('SLJ TITLE REQUIRED  ', 'ncrread   ');
  title2:=title;
  addsuffix(title2, sljinfosuffix); 
  tabout; write('READING SLJ INFO FROM ');  
  writetitle(output, title2); writeln;
  checkpresent(title2);
  textreset(textin, title2);
  title2:=title;
  addsuffix(title2, mebinsuffix); 
  tabout; write('READING OLD SLJ MES FROM ');  
  writetitle(output, title2); writeln;
  checkpresent(title2);
  rmereset(diskmein, title2);

  if not gettitle(crosstitle) then 
      warn('NO CROSSWHITE SLJ TI', 'ncrread   ');
  tabout; write('READING CROSSWHITE SLJ INFO FROM ');  
  writetitle(output, crosstitle); writeln;
  checkpresent(crosstitle);
  textreset(textin2, crosstitle);

  if not gettitle(title) then warn('SLJ OP TITLE REQUIRE', 'ncrread   ');
  title2:=title;
  addsuffix(title2, sljinfosuffix); 
  tabout; write('WRITING SLJ INFO TO ');  
  writetitle(output, title2); writeln;
  textrewrite(textout, title2);
  title2:=title;
  addsuffix(title2, mebinsuffix); 
  tabout; write('WRITING SLJ MATRIX ELEMENTS TO ');  
  writetitle(output, title2); writeln;
  rmerewrite(diskmeout, title2);  

  { read state information }
  
  getfilecard(textin,buff); tabout; write(buff); writeln;
  bcol:=1; nelectrons:=readint(buff,bcol); 
  getfilecard(textin,buff); tabout; write(buff); writeln;
  getfilecard(textin,buff);
  bcol:=1;  nslj:=readint(buff,bcol);
  if nslj>maxnslj then warn('increase maxnslj    ', 'ncrread   ');
  skipblanks(buff,bcol); skiptoblank(buff,bcol);
  lengthslj:=readint(buff,bcol);
  if lengthslj>maxsljlab then warn('increase maxsljlab  ', 'ncrread   ');
  writeln; tabout;
  write(nslj:5,' STATES');
  for i:=7 to maxsljlab do write(' ');
  write('   Y  2S  2L   X  2J WHICHSL');
  writeln; writeln;
  for wslj:=1 to nslj do
  with sljlab[wslj] do begin
    getfilecard(textin,buff); bcol:=1;
    for i:=1 to maxsljlab do sljpretty[i]:=buff[i];
    bcol:=maxsljlab+1; {skip past label}
    sljy:=readint(buff,bcol);
    sljs:=readint(buff,bcol);
    sljl:=readint(buff,bcol);
    sljx:=readint(buff,bcol);
    sljj:=readint(buff,bcol);
    sljwhichsl:=readint(buff,bcol);
    tabout; 
    write(wslj:5,' ',sljpretty,sljy:6,sljs:6,sljl:6,sljx:6,sljj:6,
         sljwhichsl:6);
    writeln;
  end;
  writeln;
  prstatus(output);

  { write state information out to file}
  write(textout,newnelectrons:3,' ELECTRONS: SLJ STATES, MES AND RMES');
  writeln(textout);
  writecreationdate(textout, progname);
  (* tabout; write(nslj:5,' STATES');
  for i:=7 to maxsljlab do write(' ');
  write('   Y  2S  2L   X  2J WHICHSL'); writeln; writeln;*)
  write(textout,nslj:4,'  STATES, ',lengthslj:4,'  CHARACTER LABELS'); 
  writeln(textout);
  for wslj:=1 to nslj do
    with sljlab[wslj] do begin
      (* tabout; 
      write(wslj:5,' ',sljpretty,sljy:6,sljs:6,sljl:6,sljx:6,sljj:6
          ,sljwhichsl:6); writeln; *)
      write(textout,sljpretty, sljy:6,sljs:6,sljl:6,sljx:6,sljj:6
        ,sljwhichsl:6); writeln(textout)
    end;
  writeln;

  { read and write the matrix elements }
  finished :=false;
  while not finished do begin
    getfilecard(textin,buff); bcol:=1;
    getstr(buff,bcol,buff1);
    for i:=1 to 10 do lab[i]:=buff1[i];
    ts:=readint(buff,bcol);
    tl:=readint(buff,bcol);
    tj:=readint(buff,bcol);
    getfilecard(textin,buff); bcol:=1; nme:=readint(buff,bcol);
    writeln(lab, ts:6, tl:6, tj:6, nme:10);
    writeln(textout, lab, ts:6, tl:6, tj:6);
    if tj<>0 then begin {only search for scalar operators}
      replacematrix:=false;
    end else begin
      replacematrix:=false; 
      textreset(textin2, crosstitle);
      if  getmatrix(lab, textin2) then replacematrix:=true;
    end; 
    if replacematrix then begin
      tabout; writeln('SUBSTITUTING NEW MATRIX ELEMENTS');
      for wme:=1 to nme do rmeread(diskmein,sljme);
      newnme:=0;
      for statei:=1 to nslj do
        for statej:=statei to nslj do
          if sljmatrix[statei,statej]<>0 then begin
            newnme:=newnme+1;
            sljme.mi:=statei; sljme.mj:=statej; 
            sljme.me:=sljmatrix[statei,statej];
            rmewrite(diskmeout,sljme);
          end;
      writeln(textout,newnme);
      tabout; writeln('NUMBER OF MATRIX ELEMENTS IS: ',newnme:1);
    end else begin
      tabout; writeln('USING OLD MATRIX ELEMENTS');
      for wme:=1 to nme do begin
        rmeread(diskmein,sljme);
        rmewrite(diskmeout,sljme);
      end;
      writeln(textout,nme);
      tabout; writeln('NUMBER OF MATRIX ELEMENTS IS: ',nme:1);
    end; 

    { special case: get T2 for f12}
    if (newnelectrons=12) and (lab='GAMMA     ') then begin
      tabout; writeln('LOOKING FOR T2');
      lab:='T2        ';
      ts:=0; tl:=0; tj:=0;     
      writeln(lab, ts:6, tl:6, tj:6, nme:10);
      writeln(textout, lab, ts:6, tl:6, tj:6);
      textreset(textin2, crosstitle);
      if  getmatrix(lab, textin2) then begin
        tabout; writeln('ADDING T2 MATRIX ELEMENTS');
        newnme:=0;
        for statei:=1 to nslj do
          for statej:=statei to nslj do
            if sljmatrix[statei,statej]<>0 then begin
              newnme:=newnme+1;
              sljme.mi:=statei; sljme.mj:=statej; 
              sljme.me:=sljmatrix[statei,statej];
              rmewrite(diskmeout,sljme);
            end;
        writeln(textout,newnme);
        tabout; writeln('NUMBER OF MATRIX ELEMENTS IS: ',newnme:1);
      end else begin
         warn('COULD NOT FIND T2   ', 'ncrread   ');
      end;
    end;

    if eof(textin) then finished:=true;
  end;

  textclose(textin);
  textclose(textin2);
  textclose(textout);
  rmeclose(diskmein);
  rmeclose(diskmeout);  

end;

{$include "racahout.zp"}  

procedure initialize;
begin
  iuu:=0; ivv:=0; juu:=0; jvv:=0;
  findmaxrange;
end;

procedure examineinstructions;
var
    j1,j2,j3,j4,j5,j6,j7,j8,j9:integer;
    s1,l1,s2,l2,ts,tl,tj:integer;
    rme,me:real;
begin
  while interpreter('-', progname) do begin
    if samestring(inbuff,incol,'SLJCALC   ',8) then sljcalculation else
       
    if samestring(inbuff,incol,'NCRREAD   ',8) then ncrread else
    if samestring(inbuff,incol,'RACAH     ',6) then racahoutput else          

    if samestring(inbuff,incol,'RANGE     ',6) then setrange else
    if samestring(inbuff,incol,'MAKECFITIP',8) then makecfitinput else
    if samestring(inbuff,incol,'SIXJ      ',5) then begin
      j1:=nextint; j2:=nextint; j3:=nextint;
      j4:=nextint; j5:=nextint; j6:=nextint;
      tabout; write('SIXJ : ',j1:3,j2:3,j3:3,j4:6,j5:3,j6:3,'  =  ');
      write(sixj(j1,j2,j3,j4,j5,j6)); writeln; writeln;
    end else
    if samestring(inbuff,incol,'NINEJ     ',6) then begin
      j1:=nextint; j2:=nextint; j3:=nextint;
      j4:=nextint; j5:=nextint; j6:=nextint;
      j7:=nextint; j8:=nextint; j9:=nextint;
      tabout; write('NINEJ: ',j1:3,j2:3,j3:3,j4:6,j5:3,j6:3,j7:6,j8:3,j9:3);
      write('  =  ',ninej(j1,j2,j3,j4,j5,j6,j7,j8,j9)); writeln; writeln;
    end else
    if samestring(inbuff,incol,'T0KK      ',5) then begin
      tl:=nextint; rme:=nextnum;
      s1:=nextint; l1:=nextint; j1:=nextint;
      s2:=nextint; l2:=nextint; j2:=nextint;
      me:=t0kkval(tl,rme,s1,l1,j1,s2,l2,j2);
      tabout; write('L REDUCED ME = ',rme); writeln;
      tabout; write('< ',s1:3,l1:3,j1:3,' ]] T  0',tl:3,tl:3);
      write(' [[ ',s2:3,l2:3,j2:3,' >  = ',me); writeln; writeln;
    end else
    if samestring(inbuff,incol,'TK0K      ',5) then begin
      ts:=nextint; rme:=nextnum;
      s1:=nextint; l1:=nextint; j1:=nextint;
      s2:=nextint; l2:=nextint; j2:=nextint;
      me:=tk0kval(ts,rme,s1,l1,j1,s2,l2,j2);
      tabout; write('S REDUCED ME = ',rme); writeln;
      tabout; write('< ',s1:3,l1:3,j1:3,' ]] T ',ts:3,0:3,ts:3);
      write(' [[ ',s2:3,l2:3,j2:3,' >  = ',me); writeln; writeln;
    end else
    if samestring(inbuff,incol,'TKK0      ',5) then begin
      ts:=nextint; rme:=nextnum;
      s1:=nextint; l1:=nextint; j1:=nextint;
      s2:=nextint; l2:=nextint; j2:=nextint;
      me:=tkk0val(ts,rme,s1,l1,j1,s2,l2,j2);
      tabout; write('SL REDUCED ME = ',rme); writeln;
      tabout; write('< ',s1:3,l1:3,j1:3,' ] T ',ts:3,ts:3,0:3);
      write(' [ ',s2:3,l2:3,j2:3,' >  = ',me); writeln; writeln;
    end else
    if samestring(inbuff,incol,'TSLJ      ',5) then begin
      ts:=nextint; tl:=nextint; tj:=nextint; rme:=nextnum;
      s1:=nextint; l1:=nextint; j1:=nextint;
      s2:=nextint; l2:=nextint; j2:=nextint;
      me:=tsljval(ts,tl,tj,rme,s1,l1,j1,s2,l2,j2);
      tabout; write('SL REDUCED ME = ',rme); writeln;
      tabout; write('< ',s1:3,l1:3,j1:3,' ]] T ',ts:3,tl:3,tj:3);
      write(' [[ ',s2:3,l2:3,j2:3,' >  = ',me); writeln; writeln;
    end else
    if samestring(inbuff,incol,'TMAGMOM   ',5) then begin
      s1:=nextint; l1:=nextint; j1:=nextint;
      s2:=nextint; l2:=nextint; j2:=nextint;
      me:=tmagmomval(s1,l1,j1,s2,l2,j2);
      tabout; write('< ',s1:3,l1:3,j1:3,' ]] L + GS');
      write(' [[ ',s2:3,l2:3,j2:3,' >  = ',me); writeln; writeln;
    end else warn('INVALID INSTRUCTION ', 'examineins');
  end;
end {examineinstructions};

procedure tidyup {tidying up procedure - not needed right now};
begin end;

procedure progmain;
begin {sljcalc}
  startup(progname);  setfactorial; initialize;
  examineinstructions;
  tidyup; finish(false);
end {sljcalc};

begin
  progmain
end.
"""
