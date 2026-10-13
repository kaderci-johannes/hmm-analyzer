void Progbar(int i, int N, string s="", int L=50)
{
	i==0 ? cout<<"\n\n\n\n" : cout;
	i++;
	std::cout<<"\033[F\033[F\033[F\033[F  Analyzing data: "<<s<<"\n    [";
	for(int j=0;j<i*L/N;j++) std::cout<<"=";
	i<N ? std::cout<<">" : cout;
	for(int j=0;j<L-i*L/N-1;j++) std::cout<<" ";
	std::cout<<std::fixed<<std::setprecision(2)<<"]  "<<100.*i/N<<"%  "<<"("<<i<<" / "<<N<<")\n\n\n";
	std::cout.flush();
	return;
}
