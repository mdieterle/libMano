
void test_tfile()
{

gSystem->Load("libMano.so");

TFile* a = new TFile("a.root");
TFile* b = new TFile("b.root");
TFile* c = new TFile("c.root");
TFile* d = new TFile("d.root");
TFile* e = new TFile("e.root");

TFile** f
f = new TFile*[3];
f[0] = c->Clone("c_0");
f[1] = d->Clone("c_1");
f[2] = e->Clone("c_2");

TMObjectCollection* m = new TMObjectCollection("m","m");
m->AddArrayObject(3,f);
//m->AddBackgroundObject(c);
m->AddSignalObject(b);
m->AddDataObject(a);

TMCuts::FitMM(m);
}
