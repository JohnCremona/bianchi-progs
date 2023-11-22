// FILE QIDLOOP.CC

#include "qidloop.h"

Qidealooper::Qidealooper(long nmin, long nmax, int both_conjugates, int sorted)
  : n(nmin), maxn(nmax), both(both_conjugates), sort(sorted)
{
  vector<Qideal> Ilist = (sort? Qideal_lists::ideals_with_norm(n, both): ideals_with_norm(n, both));
  //  cout<<" Qidealooper inserting "<<Ilist.size()<<" ideals of norm "<<n<<" into queue: "<<Ilist<<endl;
  I_norm_n.insert(I_norm_n.end(), Ilist.begin(), Ilist.end());
  advance();
}

Qideal Qidealooper::next() // must only call when not not_finished
{
  Qideal I = I_norm_n.front();
  I_norm_n.pop_front();
  advance();
  return I;
}

void Qidealooper::advance()
{
  if (! I_norm_n.empty()) return;
  while (I_norm_n.empty() && n<maxn)
    {
      n+=1;
      vector<Qideal> Ilist = (sort? Qideal_lists::ideals_with_norm(n, both): ideals_with_norm(n, both));
      //      cout<<" Qidealooper inserting "<<Ilist.size()<<" ideals of norm "<<n<<" into queue: "<<Ilist<<endl;
      I_norm_n.insert(I_norm_n.end(), Ilist.begin(), Ilist.end());
    }
  // now either n>maxn and I_norm_n is empty, or n<=maxn and it is
  // not.  In the first case, not_finished() will be false
}

vector<Qideal> Quad::class_group;                // list of ideals representing ideal classes (no structure)
vector<Qideal> Quad::class_group_2_torsion;      // list of ideals representing 2-torsion in class group
vector<Qideal> Quad::class_group_2_torsion_gens; // list of ideals generating 2-torsion in class group
vector<Qideal> Quad::class_group_2_cotorsion;      // list of ideals representing class group mod squares
vector<Qideal> Quad::class_group_2_cotorsion_gens; // list of ideals generating class group mod squares

void Quad::fill_class_group()
{
  // NB class_number may already be set correctly (if <=5 currently); if not, it will have been set to 0

  class_group.clear();
  class_group_2_torsion.clear();
  class_group_2_cotorsion.clear();
  class_group_2_torsion_gens.clear();
  class_group_2_cotorsion_gens.clear();

  class_group.push_back(Qideal());
  class_group_2_torsion.push_back(Qideal());
  class_group_2_cotorsion.push_back(Qideal());

  if(class_number==1)
    {
      class_group_2_rank = 0;
      return;
    }

  long MB = floor(2*sqrt(I2long(absdisc))/PI); // Minkowski bound
  Qidealooper loop(1, MB, 1, 1);
  while( loop.not_finished() )
    {
      Qideal I = loop.next();
      if (I.is_principal())
        continue;
      if (find_ideal_class(I, class_group) == -1)
        {
          // cout << I << " is in a new ideal class (#" << class_group.size() << ")" << endl;
          class_group.push_back(I);
          Qideal I2 = I*I;
          // cout << " -- the square is " << I2 << endl;
          if (I2.is_principal())
            {
              // cout << " -- this class is 2-torsion" <<endl;
              // I is 2-torsion
              if (find_ideal_class(I, class_group_2_torsion)==-1)
                {
                  // add a new 2-torsion generator
                  class_group_2_torsion_gens.push_back(I);
                  // compute the new coset
                  vector<Qideal> new_2_torsion;
                  for (vector<Qideal>::iterator J=class_group_2_torsion.begin(); J!=class_group_2_torsion.end(); J++)
                    new_2_torsion.push_back(I*(*J));
                  // append the new coset
                  class_group_2_torsion.insert(class_group_2_torsion.end(), new_2_torsion.begin(), new_2_torsion.end());
                }
            }
          else
            { // I has order>1, and we keep its conjugate (in the inverse class) too
              // cout << " -- this class has order > 2" <<endl;
              class_group.push_back(I.conj());
            }
        }
      if (class_group.size()==(unsigned)class_number)
        break;
    }
  if (class_number==0) // then it wasn't preset
    {
      class_number = class_group.size();
    }
  // cout << "Class number = " << class_number << " with representatives " << class_group << endl;

  // replace the 2-torsion class representatives with the equivalent
  // ideals in the main class_group list
  for (vector<Qideal>::iterator J=class_group_2_torsion.begin(); J!=class_group_2_torsion.end(); J++)
    {
      *J = class_group[find_ideal_class(*J, class_group)];
    }
  assert (class_group_2_torsion_gens.size()+1 == prime_disc_factors.size());
  class_group_2_rank = class_group_2_torsion_gens.size();

  // find ideals generating and representing Cl/Cl^2
  for (vector<Qideal>::iterator Ii=class_group.begin()+1; Ii!=class_group.end(); ++Ii)
    {
      Qideal I = *Ii;
      if (!I.has_square_class()) // then I is not a square
        {
          if (find_ideal_class_mod_squares(I, class_group_2_cotorsion)==-1)
            {
              // add a new 2-cotorsion generator
              class_group_2_cotorsion_gens.push_back(I);
              // compute the new coset
              vector<Qideal> new_2_cotorsion;
              for (vector<Qideal>::iterator Ji=class_group_2_cotorsion.begin(); Ji!=class_group_2_cotorsion.end(); ++Ji)
                {
                  Qideal J = *Ji;
                  new_2_cotorsion.push_back(class_group[find_ideal_class(I*J, class_group)]);
                }
              // append the new coset
              class_group_2_cotorsion.insert(class_group_2_cotorsion.end(), new_2_cotorsion.begin(), new_2_cotorsion.end());
            }
        }
    }
}


// END OF FILE QIDLOOP.CC
