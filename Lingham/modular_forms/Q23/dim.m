load "core.m";
load "hecke.m";
load "levels.m";
load "bianchi.m";

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
   for p in primes do
      HM cat:=[*Hecke(p,L,info[1],info[2])*];
   end for;

   if nd eq 0 then
      return "";
   elif nd eq 1 then
      time x,y:=Eigen(L,info,ns,primes,HM);
      printf "Fricke:%o\n",y;
      x;
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

//Main(50,6,6);
L := ideal<OK|[6,w]>;
Main_L(PrimeIdeals(50),L);
