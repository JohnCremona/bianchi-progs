function PrimeIdeals(nbound)

   list:=[];
   X:=[x: x in [1..nbound] | IsPrime(x)];
   for y in X do
      I:=ideal<OK|y>;
      F:=Factorization(I);
      if #F eq 1 then
         if Norm(F[1][1]) le nbound then
            list cat:=[F[1][1]];
         end if;
      else
         list cat:=[F[1][1],F[2][1]];
      end if;
   end for;

   Sort(~list,func<x,y|10*(Norm(x)-Norm(y))+cm(hinv(x))-cm(hinv(y))>);

   return list;

end function;

function Levels(nbound)

   T:=PrimeIdeals(nbound);   

   prime:=T[1];
   count:=[0..Ceiling(Log(Norm(prime),nbound))];
   store:=[prime^p:p in count];

   T:=[T[i]:i in [2..#T]];

   for prime in T do
      for s in store do
         power:=Ceiling(Log(Norm(prime),nbound)-Log(Norm(prime),Norm(s)));
         tempcount:=[1..power];
	 tempstore:=[s*prime^p:p in tempcount];
	 store cat:=tempstore;
      end for;
   end for;

   Z:=[store[i]:i in [2..#store]];

   Z:=[z:z in Z | Norm(z) le nbound];

   Sort(~Z,func<x,y|Norm(x)-Norm(y)>);

   return Z;

end function;
