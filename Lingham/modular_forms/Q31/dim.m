load "core3.m";
load "hecke3.m";
load "levels.m";
load "bianchi3.m";

E:=EllipticCurve([1,-w-1,1,w-3,1]);

function IdealClass(n,L);
   if n eq 1 then   
      if IsSquare(Norm(L)) or Norm(L) eq 31 then
         return "";
      else
         tf,g:=IsPrincipal(L);
         if g[1]*g[2] gt 0 then
            return "a";
         else
            return "b";      
         end if;
      end if;
   elif n eq 2 then
      return "a";
   elif n eq 3 then
      return "b";
   end if;
end function;

function MyNorm(I)  
   tf,sr:=IsPower(Norm(I));
   if tf eq true then
      return sr;
   else
      return Norm(I);
   end if;
end function;

function Main_L(primes,L)

   F:=Factorization(L);
   printf "Norm:%o\n",Norm(L);
   "Factorization:";
   for f in F do
      ic:=IdealClass(cm(hinv(f[1])),f[1]);
      if f[2] eq 1 then
         printf "p%o%o\n",MyNorm(f[1]),ic;
      else
         printf "p%o%o^%o\n",MyNorm(f[1]),ic,f[2];
      end if;
   end for;

   time info:=DimH(L,1);     

   fd:=Dimension(info[5]);  

   printf "Dim=%o\n",fd;          // dimension of the full homology space

   if IsPrime(L) eq true then
      ns:=info[5];
      nd:=fd;
   elif fd eq 0 then
      nd:=0;
      return "";
   elif fd eq 1 then
      nd:=1;
      ns:=info[5];
   else
      time ns:=NewS(L,info);     
      nd:=Dimension(ns);
   end if;
     
   if fd ge 1 then
      printf "NewDim=%o\n",nd;    // dimension of the newspace
   end if;

   HM:=[**];           // compile a list of all the operators in terms
                       // of matrices
//   for p in primes do
//       tp:=Hecke(p,L,info[1],info[2]);
//       HM cat:=[*tp*];
//   end for;

   if nd eq 0 then
      return "";
   elif nd eq 1 then
       for i in [1..#primes] do
	   P:=primes[i];
	   print "P=",P," of norm ",Norm(P),": ";
	   time tp:=Hecke(P,L,info[1],info[2]);
	   time x,y:=Eigen(L,info,ns,[P],[tp]);
	   if i eq 1 then printf "Fricke:%o\n",y; end if;
	   print x[1];
	   if x[1] eq TraceOfFrobenius(Reduction(E,P)) then
	       print "OK";
	   else
	       print "***WRONG***";
	   end if;
	   
       end for;
      return "";
   else                           
      time r:=Eigen(L,info,ns,primes,HM);
      for rr in r do
         printf "Fricke:%o\n",rr[2];
         rr[1];
         "";      
      end for;
      return "";
   end if;

end function;

function Main(pbound,lbound_lb,lbound_ub)

   if lbound_lb eq 0 then 
      Levels_lb:=[];
   else
      Levels_lb:=Levels(lbound_lb-1);  
   end if;

   Levels_ub:=Levels(lbound_ub);    
   Levels2:=[l:l in Levels_ub | l notin Levels_lb];
   primes:=PrimeIdeals(pbound);

   for L in Levels2 do
      Main_L(primes,L);
   end for;

   return "";

end function;

//Main(10,32,32);

// Pacetti's example 1:

L:=8*OK+4*w*OK;  // = p2a^3*p2b^2;
plist:=[3, 5, 7, 11, 13, 17, 19, 23, 29, 37, 41, 43, 47, 53, 59, 67,
    71, 73, 79, 89, 109, 127, 131, 149, 173, 193, 227, 283, 293, 349,
    379, 431, 521, 577, 607, 653, 839, 857, 1031, 1063, 1117, 1303,
    1451, 1493, 1619, 1741, 2003, 2153, 2333, 2707, 2767, 2963, 3119,
    3373, 3767];
Plist:=&cat([[f[1]: f in Factorization(p*OK)] : p in plist]);

SetLogFile("pacetti-31-b.out2");
L:=Conductor(E);
plist:=[3, 7, 11, 13, 17, 19, 23, 29, 47, 59, 67, 71, 89, 97, 101, 103,
    107, 109, 149, 157, 163, 191, 193, 211, 293, 311, 317, 359, 443,
    577, 607, 617, 653, 691, 701];
Plist:=&cat([[f[1]: f in Factorization(p*OK)] : p in plist]);
Main_L(Plist,L);
UnsetLogFile();
