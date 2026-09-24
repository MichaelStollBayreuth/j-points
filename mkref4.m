// Makes testbase4, the reference for test4.sh, with Magma: for every curve of
// testcurves2 the points up to the height bound 30, found by brute force
// over every coprime triple (a, b, c) of first coordinates -- Magma's Kummer
// package solves the Kummer equation for the fourth coordinate
// (Points(K, [a, b, c])) and lifts the points to the Jacobian (Points(J, P)),
// which is independent of j-points; Magma's own search Points(J : Bound) is
// a port of j-points and is not used.  Two lists per curve, as j-points
// prints them (coprime integers, the first nonzero one positive, the origin
// left out, sorted as "LC_ALL=C sort" sorts them): with the bound on all four
// coordinates, and with -a, the bound on the triple only.
// Run with "magma -b mkref4.m > testbase4"; about five minutes.
H := 30;
P<x> := PolynomialRing(Rationals());
function normalise(s)
  d := LCM([Denominator(c) : c in s]);
  s := [Integers() | c*d : c in s];
  g := GCD(s);
  s := [c div g : c in s];
  i := 1; while s[i] eq 0 do i +:= 1; end while;
  if s[i] lt 0 then s := [-c : c in s]; end if;
  return s;
end function;
for line in Split(Read("testcurves2"), "\n") do
  if #line eq 0 or line[1] eq "#" then continue; end if;
  fields := Split(line, ";");
  coeffs := [StringToInteger(c) : c in Split(fields[1], " ")];
  f := P!coeffs;
  J := Jacobian(HyperellipticCurve(f));
  K := KummerSurface(J);
  plain := {}; all := {};
  for a in [0..H] do for b in [-H..H] do for c in [-H..H] do
    if a eq 0 and (b lt 0 or (b eq 0 and c le 0)) then continue; end if;
    if GCD([a, b, c]) ne 1 then continue; end if;
    for pt in Points(K, [a, b, c]) do
      if #Points(J, pt) eq 0 then continue; end if;
      s := normalise(Eltseq(pt));
      str := Sprintf("(%o, %o, %o, %o)", s[1], s[2], s[3], s[4]);
      Include(~all, str);
      if Max([Abs(c) : c in s]) le H then Include(~plain, str); end if;
    end for;
  end for; end for; end for;
  printf "# j-points '%o' %o\n", fields[1], H;
  for s in Sort(Setseq(plain)) do print s; end for;
  printf "# j-points '%o' %o -a\n", fields[1], H;
  for s in Sort(Setseq(all)) do print s; end for;
end for;
quit;
