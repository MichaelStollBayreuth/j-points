// Makes testbase2, the reference for test2.sh, with Magma: for every line of
// testcurves2 the points on the Jacobian up to the height bound, found by
// Magma's own search (Points(J : Bound := h)), written as the Kummer
// coordinates j-points prints (coprime integers, the first nonzero one
// positive, the origin left out), sorted as "LC_ALL=C sort" sorts them.
// Run with "magma -b mkref2.m > testbase2".
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
  h := StringToInteger(fields[2]);
  printf "# j-points '%o' %o\n", fields[1], h;
  f := P!coeffs;
  if not IsSquarefree(f) then print "NOT SQUAREFREE"; continue; end if;
  J := Jacobian(HyperellipticCurve(f));
  K := KummerSurface(J);
  strs := {};
  for pt in Points(J : Bound := h) do
    s := normalise(Eltseq(K!pt));
    if s[1] eq 0 and s[2] eq 0 and s[3] eq 0 then continue; end if;
    Include(~strs, Sprintf("(%o, %o, %o, %o)", s[1], s[2], s[3], s[4]));
  end for;
  for s in Sort(Setseq(strs)) do print s; end for;
end for;
quit;
