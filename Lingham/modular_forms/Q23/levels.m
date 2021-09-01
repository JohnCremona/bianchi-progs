// Function to order prime ideals: first by norm, or if the norms are
// equal then by the ideal class.  This relies on cm(hinv(x)), for x
// an ideal, returning a number.

prime_ideal_sort_func := func<x,y|
Norm(x) eq Norm(y) select cm(hinv(x))-cm(hinv(y))
                   else (Norm(x)-Norm(y))>;

// Function to return a sorted list of prime ideals of norm up to nbound:
	       
function PrimeIdeals(nbound)

   plist:=[P : P in &cat([[f[1] : f in Factorization(p*OK)] : p in [2..nbound] | IsPrime(p)]) | Norm(P) le nbound]; 

   Sort(~plist,prime_ideal_sort_func);

   return plist;

end function;

// Recursive function which returns a partially sorted list of all
// integral ideals of norm up to nbound with support in plist
// (or with unrestricted support by default).
// The partial sorting is by norm, but with no tie-break in case of
// repeated norm.

function Levels(nbound :  plist:=[])
    if nbound eq 0 then return []; end if;
    if nbound eq 1 then return [1*OK]; end if;
    if plist eq [] then plist:=PrimeIdeals(nbound); end if;
    prime:=plist[1];
    Exclude(~plist,prime);
    count:=[0..Floor(Log(Norm(prime),nbound))];
    store:=[prime^p:p in count];
    if #plist eq 0 then return store; end if;
    list := &cat([[A*B : B in $$(Floor(nbound/Norm(A)) : plist:=plist)] : A in store]);
    list:=[A: A in list | Norm(A) gt 1 and Norm(A) le nbound];
    Sort(~list,func<x,y|Norm(x)-Norm(y)>);
    return list;
end function;
