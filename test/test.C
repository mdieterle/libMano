
void test()
{

gSystem->Load("libMano.so");



TH1D* h1 = new TH1D("h_0","h_0",1000,0,1000);
TH1D* h2 = new TH1D("h_1","h_1",1000,0,2000);
TH1D* h3 = new TH1D("h_2","h_2",1000,0,3000);
TH1D* h4 = new TH1D("h_3","h_3",1000,0,5000);

TH1D** h;
h = new TH1D*[3];
h[0] = h1;
h[1] = h2;
h[2] = h3;
TMObjectCollection* m = new TMObjectCollection("m","m");
m->AddArrayObject(3,h);
m->AddBackgroundObject(h4);

if (TMTools::IsObjectType(h[0],"TH1D")) cout << "good" << endl;
else gSystem->Exit(0);

Int_t n = m->GetArraySize();
TH1D** b;
TH1D** e;

b = new TH1D*[n];

b = (TH1D**)m->GetArrayObject();

e = new TH1D*[0];
e[0] = (TH1D*)m->GetBackgroundObject();

TH1D** d;
d = new TH1D*[0];
d[0] = new TH1D("d","d",1000,0,4000);

TCanvas* c = new TCanvas("c","c",800,800);
c->Divide(3,3);

c->cd(1);
b[0]->Draw();
c->cd(2);
b[1]->Draw();
c->cd(3);
b[2]->Draw();
c->cd(4);
d[0]->Draw();
c->cd(5);
e[0]->Draw();

}
