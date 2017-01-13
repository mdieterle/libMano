void keylist()
{
   TFile f("/usr/users/dieterle/Macros/MissMass/gaga.root");
   TIter next(f.GetListOfKeys());
   TKey *key;
   while ((key=(TKey*)next())) {
//      printf("key: %s points to an object of class: %s at %d\n",
//      key->GetName(),
//      key->GetClassName(),key->GetSeekKey());
      if (strcmp(key->GetClassName(),"TOHistoCollection") != 0)
         printf("No object of class TOHistoCollection found\n");
      else
      {
         printf("Found object of class TOHistoCollection of name %s\n", key->GetName());
         ((TOHistoCollection*)f->Get(key->GetName())).Print();
      }
   }
}
