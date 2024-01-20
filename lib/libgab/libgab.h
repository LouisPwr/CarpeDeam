#ifndef libgab_h
#define libgab_h


#include <stdint.h>
#include <bitset>
#include <cmath> 
#include <limits> 
#include <string> 
#include <vector> 
#include <stdio.h>
#include <stdlib.h>
#include <sstream> 
#include <iostream>
#include <fstream>
#include <ctime>
#include <map> 
#include <iomanip>
#include <memory>
#include <algorithm> 
#include <functional> 
#include <sys/time.h> //for srand
#include <iterator>


//For directory/file operations
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h> 
#include <dirent.h>
#include <gzstream.h>
#include <errno.h>
#include <sys/resource.h>
#include <pwd.h>


using namespace std;


inline char upper(const char c){
    return char(toupper(c));
}


//Check if it is either A,C,G,T,N
inline bool isValidDNA(const char c){
    char _c= upper(c);

    if(_c ==    'A')
	return true;
    if(_c ==    'C')
	return true;
    if(_c ==    'G')
	return true;
    if(_c ==    'T')
	return true;
    if(_c ==    'N')
	return true;
    return false;
}



//Check if it is either on of the ambiguous UIPAC code
inline bool isAmbiguousUIPAC(const char c){
    char _c= upper(c);

    if(_c ==    'R')
	return true;
    if(_c ==    'Y')
	return true;
    if(_c ==    'S')
	return true;
    if(_c ==    'W')
	return true;
    if(_c ==    'K')
	return true;
    if(_c ==    'M')
	return true;
    if(_c ==    'B')
	return true;
    if(_c ==    'D')
	return true;
    if(_c ==    'H')
	return true;
    if(_c ==    'V')
	return true;

    return false;
}




//Check if it is either A,C,G,T
inline bool isResolvedDNA(const char c){
    char _c= upper(c);

    if(_c ==    'A')
	return true;
    if(_c ==    'C')
	return true;
    if(_c ==    'G')
	return true;
    if(_c ==    'T')
	return true;
    return false;
}

inline int base2int(const char c){
    char _c= upper(c);

    if(_c ==    'N')
	return 0;
    if(_c ==    'A')
	return 1;
    if(_c ==    'C')
	return 2;
    if(_c ==    'G')
	return 3;
    if(_c ==    'T')
	return 4;
    cerr<<"libgab.h base2int() Invalid base "<<c<<endl;
    exit(1);
}

inline int baseResolved2int(const char c){
    char _c= upper(c);
    if(_c == 'A')
	return 0;
    if(_c == 'C')
	return 1;
    if(_c == 'G')
	return 2;
    if(_c == 'T')
	return 3;
    cerr<<"libgab.h baseResolved2int() Invalid base "<<c<<endl;
    exit(1);
}

//Returns an index for every 2mer of different base A,C,G,T
inline int dimer2index(const char c1,const char c2){
    char _c1= upper(c1);
    char _c2= upper(c2);


    if(_c1     ==    'A'){

	if(_c2 ==    'C')
	    return 0;
	if(_c2 ==    'G')
	    return 1;
	if(_c2 ==    'T')
	    return 2;

	cerr<<"Libgab.h:1 dimer2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }


    if(_c1     ==    'C'){

	if(_c2 ==    'A')
	    return 3;
	if(_c2 ==    'G')
	    return 4;
	if(_c2 ==    'T')
	    return 5;

	cerr<<"Libgab.h:2 dimer2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }


    if(_c1     ==    'G'){

	if(_c2 ==    'A')
	    return 6;
	if(_c2 ==    'C')
	    return 7;
	if(_c2 ==    'T')
	    return 8;

	cerr<<"Libgab.h:3 dimer2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }



    if(_c1     ==    'T'){

	if(_c2 ==    'A')
	    return 9;
	if(_c2 ==    'C')
	    return 10;
	if(_c2 ==    'G')
	    return 11;

	cerr<<"Libgab.h:4 dimer2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }



    cerr<<"Libgab.h:5 dimer2index invalid dimer "<<c1<<" "<<c2<<endl;
    exit(1);
}



//Returns an index for every 2mer of different base A=0,C=1,G=2,T=3
inline int dimer2indexInt(const int c1,const int c2){
		
    if(c1     ==    0){

	if(c2 ==    1)
	    return 0;
	if(c2 ==    2)
	    return 1;
	if(c2 ==    3)
	    return 2;

	cerr<<"Libgab.h:1 dimer2indexInt invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }


    if(c1     ==    1){

	if(c2 ==    0)
	    return 3;
	if(c2 ==    2)
	    return 4;
	if(c2 ==    3)
	    return 5;

	cerr<<"Libgab.h:2 dimer2indexInt invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }


    if(c1     ==    2){

	if(c2 ==    0)
	    return 6;
	if(c2 ==    1)
	    return 7;
	if(c2 ==    3)
	    return 8;

	cerr<<"Libgab.h:3 dimer2indexInt invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }



    if(c1     ==    3){//3

	if(c2 ==    0)
	    return 9;
	if(c2 ==    1)
	    return 10;
	if(c2 ==    2)
	    return 11;

	cerr<<"Libgab.h:4 dimer2indexInt invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }



    cerr<<"Libgab.h:5 dimer2indexInt invalid dimer "<<c1<<" "<<c2<<endl;
    exit(1);
}




//Returns an index for every 2mer of different
inline int twoBases2index(const char c1,const char c2){
    char _c1= upper(c1);
    char _c2= upper(c2);

    if(_c1     ==    'A'){

	if(_c2 ==    'A')
	    return 0;
	if(_c2 ==    'C')
	    return 1;
	if(_c2 ==    'G')
	    return 2;
	if(_c2 ==    'T')
	    return 3;

	cerr<<"Libgab.h:1 twoBases2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }


    if(_c1     ==    'C'){


	if(_c2 ==    'A')
	    return 4;
	if(_c2 ==    'C')
	    return 5;
	if(_c2 ==    'G')
	    return 6;
	if(_c2 ==    'T')
	    return 7;

	cerr<<"Libgab.h:2 twoBases2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }


    if(_c1     ==    'G'){

	if(_c2 ==    'A')
	    return 8;
	if(_c2 ==    'C')
	    return 9;
	if(_c2 ==    'G')
	    return 10;
	if(_c2 ==    'T')
	    return 11;


	cerr<<"Libgab.h:3 twoBases2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }



    if(_c1     ==    'T'){

	if(_c2 ==    'A')
	    return 12;
	if(_c2 ==    'C')
	    return 13;
	if(_c2 ==    'G')
	    return 14;
	if(_c2 ==    'T')
	    return 15;

	cerr<<"Libgab.h:4 twoBases2index invalid dimer "<<c1<<" "<<c2<<endl;
	exit(1);
    }



    cerr<<"Libgab.h:5 twoBases2index invalid dimer "<<c1<<" "<<c2<<endl;
    exit(1);
}




inline char complement(const char c){
    if(c ==    'A')
	return 'T';

    if(c ==    'C')
	return 'G';

    if(c ==    'G')
	return 'C';

    if(c ==    'T')
	return 'A';



    if(c ==    'a')
	return 't';

    if(c ==    'c')
	return 'g';

    if(c ==    'g')
	return 'c';

    if(c ==    't')
	return 'a';



    if(c ==    'N')
	return 'N';

    cerr<<"Libgab.h: complement: Invalid base pair="<<c<<endl;
    exit(1);
}

//A=0,C=1,G=2,T=3
inline int complementInt(const int c){
    if( c>=0 && c<=3 )
	return (3-c);

    /* if(c ==     0) */
    /* 	return  3; */

    /* if(c ==     1) */
    /* 	return  2; */

    /* if(c ==     2) */
    /* 	return  1; */

    /* if(c ==     3) */
    /* 	return  0; */


    cerr<<"Libgab.h: complementInt: Invalid base pair="<<c<<endl;
    exit(1);
}



inline string reverseComplement(const string & inputString){
    string toReturn="";
    if(inputString.size() >0 )	
	for(int i=(inputString.size()-1);i>=0;i--){
	    toReturn+=complement( inputString[i] );
	}

    return toReturn;
}

inline bool starsWith(const string & header,const string & tocheck){
    return (header.substr(0, min(header.length(),tocheck.length()) ) == tocheck);
}



	
inline string boolStringify(const bool b){
    return b ? "true" : "false";
}


inline string booleanAsString(bool toprint){
    if(toprint)
	return string("currently turned on/used");
    else
	return string("not on/not used");
}


inline bool validOneBP(string & toCheck){
    if(toCheck.length() != 1)
	return false;
    char mychar=toupper(toCheck[0]);

    if(mychar == 'A')
	return true;
    if(mychar == 'C')
	return true;
    if(mychar == 'G')
	return true;
    if(mychar == 'T')
	return true;

    return false;
}

inline bool validAltBP(const string & toCheck){
    if(toCheck.length() != 1)
	return false;
    char mychar=toupper(toCheck[0]);

    if(mychar == '.')
	return true;
    if(mychar == 'A')
	return true;
    if(mychar == 'C')
	return true;
    if(mychar == 'G')
	return true;
    if(mychar == 'T')
	return true;

    return false;
}



/* inline string returnFirstToken(string & toparse,string & delim){ */
/*     size_t found; */
/*     string toreturn; */
/*     found=toparse.find(delim); */
/*     if (found!=string::npos){ */
/* 	toreturn=toparse.substr(0,found); */
/* 	toparse.erase(0,found+delim.length()); */
/* 	return toreturn; */
/*     } */
/*     toreturn=toparse; */
/*     toparse.erase(0,toparse.length()); */
/*     return toreturn; */
/* } */


/* inline vector<string> allTokens(string  & toparse,string & delim){ */
/*     vector<string> toReturn; */
   
/*     string token=returnFirstToken(toparse,delim); */
/*     while(token.length()!=0){ */
/* 	toReturn.push_back(token); */
/* 	token=returnFirstToken(toparse,delim); */
/*     } */
		
/*     return toReturn; */
/* } */



inline vector<string> allTokens(const string  & toparse,const char  delim){
    vector<string> toReturn;
    size_t lastfound=-1;
    while(true){
	size_t found = toparse.find(delim,lastfound+1);
	if (found!=string::npos){
	    toReturn.push_back(toparse.substr(lastfound+1,found-lastfound-1));
	    lastfound=found;
	}else{
	    toReturn.push_back(toparse.substr(lastfound+1));
	    break;
	}
    }
	       
    return toReturn;
}


inline vector<string> allTokensWhiteSpaces(const string  & toparse){

    istringstream iss(toparse);
    //solution taken from https://stackoverflow.com/questions/236129/the-most-elegant-way-to-iterate-the-words-of-a-string
    vector<string> toReturn; 
 
    copy(istream_iterator<string>(iss),
	 istream_iterator<string>(),
	 back_inserter(toReturn));

	       
    return toReturn;
}


/* inline vector<string> allTokens(const string  & toparse,const string & delim){ */
/*     if(delim.size() == 0){ */
/* 	cerr<<"Libgab.h: allTokens: delim must have at least one char"<<endl; */
/* 	exit(1); */
/*     } */
/*     /\* if(delim.size() == 1) *\/ */
/*     /\* 	return allTokens(toparse,delim[0]); *\/ */

/*     cerr<<"Libgab.h: allTokens: to implement"<<endl; */
/*     exit(1); */

/*     /\* vector<string> toReturn; *\/ */

 
		
/*     /\* return toReturn; *\/ */
/* } */


inline bool isInsert(const string & toCheck){
    vector<string> allt =allTokens(toCheck,',');
    for(unsigned int i=0;i<allt.size();i++){
	if(allt[i].length() > 1)
	    return true;
    }
    return false;
}




inline map<string, string> * info2map(const string & info){
    //auto_ptr< map<string, string> > toReturn  (new map<string, string>() );
    map<string, string> * toReturn = new map<string, string>() ;
    vector<string> fields=allTokens(info,';');
    for(unsigned int i=0;i<fields.size();i++){
	size_t found = fields[i].find("=");
	if (found!=string::npos){
	    (*toReturn)[fields[i].substr(0,found)]=fields[i].substr(found+1);
	}else{
	    (*toReturn)[fields[i]]="";
	}
    }
    return toReturn;
}




inline unsigned int string2uint(const string & s){
    return (unsigned int)strtoul(s.c_str(),NULL,0);
}



template <typename T>
string stringify(const T i){
    stringstream s;
    s << i;
    return s.str();
}
	
template <typename T>
T destringify( const string& s ){
    istringstream i(s);
    T x;
    if (!(i >> x)){
	cerr<<"Libgab.cpp: destringify() Unable to convert string=\""<<s<<"\""<<endl;
	exit(1);
    }
    return x;
} 



inline bool isInt( const string& s ){
    istringstream i(s);
    uint64_t    x;
    double      y;

    if( !(i >> x) ){
	return false;
    }
    
    //check if double
    if( i >> y ){
	if(x!=uint64_t(y))
	    return false;	
    }

    return true;
} 

inline bool isPositiveInt( const string& s ){
    istringstream i(s);
    int    x;
    double y;

    if( !(i >> x) ){
	return false;
    }
    
    //check if double
    if( i >> y ){
	if(x!=int(y))
	    return false;	
    }

    if(x<0)
	return false;

    return true;
} 
	

template <typename T>
string var2binary(const T i){
    stringstream s;
    std::bitset<8*sizeof(i)> x(i);
    s << x;
    return s.str();
}
	

inline bool isDirectory(const string & dir){
    struct stat st;
    if(stat(dir.c_str(),&st) == 0)
	if( st.st_mode & S_IFDIR )
	    return true;
    
    return false;
}


inline bool isFile(const string & dir){
    struct stat st;
    if(stat(dir.c_str(),&st) == 0)
	if( st.st_mode & S_IFREG )
	    return true;
    
    return false;
}


inline vector<string>  getdir (const string & dir){
    vector<string> toReturn;
    DIR *dp;
    struct dirent *dirp;
    if(  !( dp  = opendir(dir.c_str())  ) ){
	cerr <<"Cannot open directory = "<<dir<<endl;
	exit(1);   
    }

    while(  (dirp = readdir(dp)) ){
        toReturn.push_back(string(dirp->d_name));
    }
    closedir(dp);


    return toReturn;
}


inline int editDistTwoStrings (string const & str1, string const & str2){
    if (str1.length() != str2.length()){
	cerr <<"editDistTwoStrings(): Cannot count edit distance for two sequences that have different lengths"<<endl;
	exit(1);   	
    }

    int editDist=0;

    for(unsigned int i=0;i<str1.length();i++){
	if ( str1[i] != str2[i] ){
	    editDist++;
	}
    }
    
    return editDist;
}

inline bool strEndsWith (string const &stringToLookIn, string const &suffix){
    if (stringToLookIn.length() >= suffix.length()) 
        return (stringToLookIn.compare (stringToLookIn.length() - suffix.length(), suffix.length(), suffix) == 0 );
    return false;
}


inline bool strBeginsWith (string const &stringToLookIn, string const &prefix){
    if (stringToLookIn.length() >= prefix.length()) 
        return (stringToLookIn.compare (0, prefix.length(), prefix) == 0 );
    return false;
}

template <typename T>
inline string arrayToString(const T toPrint[] ,const int size,const string separator=","){
    if(size == 0){
    	return "";
    }
    string toReturn="";
    for(int i=0;i<(size-1);i++){
    	toReturn+=(stringify(toPrint[i])+separator);
    }
    toReturn+=(stringify(toPrint[ size -1 ]));
    return toReturn;
}

template <typename T>
inline string vectorToString(const vector<T> & toPrint,const string separator=","){
    if(toPrint.size() == 0){
	return "";
    }
    string toReturn="";
    for(int i=0;i<(int(toPrint.size())-1);i++){
	toReturn+=(stringify(toPrint[i])+separator);
    }
    toReturn+=(stringify(toPrint[ toPrint.size() -1 ]));
    return toReturn;
}


template <typename T>
inline string iteratorToString(const T & toPrint,const char * separator=","){
    typename T::const_iterator it;
    stringstream s;
    if(!toPrint.empty()){
	typename T::const_iterator itEnd= toPrint.end();
	--itEnd;
	for(it=toPrint.begin();it!=itEnd;it++){
	    s<<*it<<separator;
	}
	s<<*itEnd;
    }else{
	s<<"";
    }    
    return s.str();
}


/* template <typename T> */
/* inline string iteratorToString(const T & toPrint,const char * separator=","){ */
/*     stringstream s; */
/*     ostream_iterator<typename T::value_type> out_it (s,separator); */

/*     if(toPrint.size() > 1){ */
/* 	copy ( toPrint.begin(), toPrint.end()-1, out_it ); */
/* 	s<<toPrint.back(); */
/*     }else{ */
/* 	if(toPrint.size() == 1){ */
/* 	    s<<toPrint.back(); */
/* 	}else{ */
/* 	    s<<""; */
/* 	} */
/*     } */
/*     return  s.str(); */
/* } */


inline string zeroPad(const int numberToOutput,const int amountOfZeros){
    stringstream s;
    s  << setfill('0') << setw(amountOfZeros) <<numberToOutput;
    return s.str();
}



//A=0,C=1,G=2,T=3
inline int randomBPExceptIntTS(const int c){
    if(c ==    0)
	return 2;
    if(c ==    1)
	return 3;
    if(c ==    2)
	return 0;
    if(c ==    3)
	return 1;

	
    cerr<<"Libgab.h randomBPExceptIntTS wrong input: "<<c<<endl;
    exit(1);	
}

//! To trim white spaces
/*!
 *
 * Trims white spaces at the beginning and end of a string

  \param str : Address of string to trim
*/
static inline void trimWhiteSpacesBothEnds(string * str) { 
    str->erase(str->begin(), find_if(str->begin(), str->end(), not1(ptr_fun<int, int>(isspace))));
    str->erase(find_if(str->rbegin(), str->rend(), not1(std::ptr_fun<int, int>(isspace))).base(), str->end());// 


}


//! Returns an integer between 0 and 15  for base {A,C,G,T}x{A,C,G,T}
/*!
 * This function returns an integer according to the following table:
 * 
 * AA  0 
 * AC  1
 * AG  2
 * AT  3
 * CA  4
 * CC  5
 * CG  6
 * CT  7
 * GA  8
 * GC  9
 * GG 10
 * GT 11
 * TA 12
 * TC 13
 * TG 14
 * TT 15
*/
inline int allelePair2Int(char bp1,char bp2){
    char _bp1= upper(bp1);
    char _bp2= upper(bp2);

    if(_bp1 == 'A'){
	if(_bp2 == 'A'){
	    return 0;
	}
	if(_bp2 == 'C'){
	    return 1;
	}
	if(_bp2 == 'G'){
	    return 2;
	}
	if(_bp2 == 'T'){
	    return 3;
	}
    }

    if(_bp1 == 'C'){
	if(_bp2 == 'C'){
	    return 5;
	}
	if(_bp2 == 'T'){
	    return 7;
	}
	if(_bp2 == 'A'){
	    return 4;
	}
	if(_bp2 == 'G'){
	    return 6;
	}
    }


    if(_bp1 == 'G'){
	if(_bp2 == 'G'){
	    return 10;
	}
	if(_bp2 == 'A'){
	    return 8;
	}
	if(_bp2 == 'C'){
	    return 9;
	}
	if(_bp2 == 'T'){
	    return 11;
	}
    }


    if(_bp1 == 'T'){
	if(_bp2 == 'T'){
	    return 15;
	}
	if(_bp2 == 'G'){
	    return 14;
	}
	if(_bp2 == 'A'){
	    return 12;
	}
	if(_bp2 == 'C'){
	    return 13;
	}
    }


    cerr<<"Libgab.h allelePair2Int invalid 2 bp: "<<bp1<<" and "<<bp2<<endl;
    exit(1);

}



inline int isPotentialTransition(const char bp1,const char bp2){
    char _bp1= upper(bp1);
    char _bp2= upper(bp2);

    if(_bp1 == 'A'){
	if(_bp2 == 'G'){
	    return true;
	}
	return false;
    }

    if(_bp1 == 'C'){
	if(_bp2 == 'T'){
	    return true;
	}
	return false;
    }


    if(_bp1 == 'G'){
	if(_bp2 == 'A'){
	    return true;
	}
	return false;
    }


    if(_bp1 == 'T'){
	if(_bp2 == 'C'){
	    return true;
	}
	return false;
    }


    cerr<<"Libgab.h isPotentialTransition invalid 2 bp: "<<bp1<<" and "<<bp2<<endl;
    exit(1);

}


inline string getDateString(){
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    return  
	""+
	stringify(now->tm_year + 1900)+"-"+
	zeroPad(  (now->tm_mon + 1) ,2)+"-"+
	zeroPad(  (now->tm_mday + 1) ,2);

}

inline string getTimeString(){
    time_t t = time(0);   // get time now
    struct tm * now = localtime( & t );
    return  
	""+
	zeroPad(now->tm_hour,2)+":"+
	zeroPad(now->tm_min ,2)+":"+
	zeroPad(now->tm_sec ,2);
}





/* static inline void trimWhiteSpacesBothEnds(string * str) { */
/*     stringstream trimVar; */
/*     trimVar << *str; */
/*     str->clear(); */
/*     trimVar >> *str;    */
/* } */



inline double correlation(const vector<double>& x, const vector<double>& y){
    if(x.size() == 0){
	cerr<<"Libgab.h correlation() ERROR: the size of x is zero"<<endl;
	exit(1);

    }

    if(x.size() != y.size()){
	cerr<<"Libgab.h correlation() ERROR: the size of x is not the same as the size of y"<<endl;
	exit(1);
    }

    double n  = double(x.size());
    double ex(0), ey(0),xt(0), yt(0),sxx(0),syy(0),sxy(0) ;

    for (unsigned int i = 0; i < n; i++) { 
	ex += x[i];
	ey += y[i];
    }
    /* cout<<n<<endl; */
    ex = ex/n;
    ey = ey/n;
    /* cout<<"ex " <<ex<<endl; */
    /* cout<<"ey " <<ey<<endl; */

    for (unsigned int i = 0; i < n; i++) { 
	xt = x[i] - ex;
	yt = y[i] - ey;
	/* cout<<"xt "<<xt<<endl; */
	/* cout<<"yt "<<yt<<endl; */

	sxx += xt * xt;
	syy += yt * yt;
	sxy += xt * yt;
    }

    /* cout<<"1 "<<sxx<<endl; */
    /* cout<<"2 "<<sxy<<endl; */
    /* cout<<"3 "<<syy<<endl; */

    return (sxy/(sqrt(sxx*syy)+ numeric_limits< double >::min() ));
}


inline string returnGitHubVersion(const string  programName,const string  suffixToAdd){
    //getting github version
    string directoryProgram;
    string commandPath=string(programName);
    size_t posSlash=commandPath.find_last_of("/");
    if(posSlash == string::npos){
        directoryProgram="";
    }else{
        directoryProgram=commandPath.substr(0,posSlash);
    }
    string gitFileLog=directoryProgram+"/"+suffixToAdd+"/.git/logs/HEAD";
    string gitVersion="NA";
    if(isFile(gitFileLog)){
        ifstream myFile;
        string line;
        myFile.open(gitFileLog.c_str(), ios::in);

        if (myFile.is_open()){
            while ( getline (myFile,line)){
                vector<string> vs=allTokens(line,' ');
                gitVersion=vs[1];
            }
            myFile.close();

        }else{
            cerr << "Unable to open github file "<<gitFileLog<<endl;
            return "NA";
        }
    }
    return gitVersion;
}


template <typename T>
inline vector<T>  vectorDist(const vector<T> & toEvaluate){
    vector<T> toReturn;
    T m2;
    T m1;
    T m;
    if(toEvaluate.size() <=1){
	return toReturn;
    }
    //have at least 2
    m2=toEvaluate[0];
    m1=toEvaluate[1];
    if(m2>m1){
	cerr<<"libgab.h vectorDist() vector is unsorted"<<endl;
	exit(1);
    }
    toReturn.push_back( (m1-m2)  );
         
    for(unsigned int i=2;i<toEvaluate.size();i++){
	m=toEvaluate[i];
	if(m1>m){
	    cerr<<"libgab.h vectorDist() vector is unsorted"<<endl;
	    exit(1);
	}

	if( (m1-m2) < (m-m1) ){ //if m1 
	    toReturn.push_back( (m1-m2)  );
	}else{
	    toReturn.push_back( (m-m1)  );
	}
	
	m2=m1;
	m1=m;
    }

    toReturn.push_back( (m1-m2)  );    

    return toReturn;
}


/* 
inline char * getRealpath(const char *arg){
    int pathmaxtouse=1024;//under MacOS, PATH_MAX is defined to be 1024 in /usr/include/sys/syslimits.h
    
#ifdef PATH_MAX
    pathmaxtouse=PATH_MAX;//if the max # of bytes is a path is defined
#endif
    
    char actualpath [pathmaxtouse+1];
    char * returnRealpath = realpath(arg, actualpath);
    return returnRealpath;    
}

inline char * getENV(const char *arg){
    string envpath = string(getenv("PATH"));
    vector<string> token=allTokens( envpath,':');
    for(unsigned int i=0;i<token.size();i++){
	string t=token[i]+"/"+string(arg);
	if( isFile( t )){
	    return getRealpath(t.c_str() );

	}
    }
    return NULL;
}



inline string getCWD(const char *arg){
    //string tm=string(arg);
    char * returnRealpath = getRealpath(arg);
    if(returnRealpath == NULL){
	returnRealpath=getENV(arg);
	if(returnRealpath == NULL){
	    cerr<<"libgab.h getCWD failed on  "<<*arg<<endl;
	    exit(1);
	}
    }
   
    vector<string> token=allTokens( string(returnRealpath),'/');
    token.pop_back();

    if(strEndsWith(vectorToString(token,"/"),"/")){
	return vectorToString(token,"/");
    }else{
	return vectorToString(token,"/")+"/";
    }
} */

inline string getHomeDir(){
    struct passwd *pw = getpwuid(getuid());
    const char *homedir = pw->pw_dir;
    return string(homedir);
}

/* inline string getFullPath(const string & st){
    vector<string> token=allTokens( st,'/');
    //sub the ~ for the home dir
    for(unsigned int i=0;i<token.size();i++)  
	if(token[i] == "~")  
	    token[i]=getHomeDir();
    string stT=vectorToString(token,"/");
    int pathmaxtouse=1024;//under MacOS, PATH_MAX is defined to be 1024 in /usr/include/sys/syslimits.h
    
#ifdef PATH_MAX
    pathmaxtouse=PATH_MAX;//if the max # of bytes is a path is defined
#endif
    
    char actualpath [pathmaxtouse+1];
    char * returnRealpath = realpath(stT.c_str(), actualpath);
    
    if(returnRealpath == NULL){
	cerr<<"libgab.h getFullPath failed on "<<st<<endl;
	exit(1);
    }

    if(isDirectory(actualpath)){
	token=allTokens( string(actualpath),'/');
    
	if(strEndsWith(vectorToString(token,"/"),"/")){
	    return vectorToString(token,"/");
	}else{
	    return vectorToString(token,"/")+"/";
	}
    }else
	return actualpath;
} */


inline pair<double,double> computeMeanSTDDEV(const vector<double> & v){
    double sum = 0.0;
    for(unsigned int i=0;i<v.size();i++)
	sum += v[i];

    double m =  sum / double(v.size());

    double accum = 0.0;

    for(unsigned int i=0;i<v.size();i++)
	accum += ( (v[i] - m) * (v[i] - m) );


    double stdev = sqrt( accum / double(v.size()-1) );
    return make_pair(m,stdev);
}


inline uint64_t nChoosek( uint64_t n, uint64_t k ){
    if (k > n) return 0;
    if (k * 2 > n) k = n-k;
    if (k == 0) return 1;
    
    uint64_t result = n-k+1;
    for(uint64_t i = 2; i <= k; ++i ) {
	result *= (n-k+i);
	result /= i;
    }

    return result;
}

//to convert a DNA string (A,C,G,T) of less than 32 characters to an unsigned 64 bits integer
inline uint64_t seq2uint64(string & s){
    uint64_t toreturn = 0;
    for(unsigned int i=0;i<s.size();i++){
	toreturn =	toreturn<<2 |( (s[i] == 'A')?0:((s[i] == 'C')?1:((s[i] == 'G')?2:3 ) ));
    }

    return toreturn;
}

/* template <typename T> */
/* template <typename T2> */
template<class T, class T2>
inline vector<T>  allKeysMap(map<T,T2> & m){
    vector<T> v;
    for(typename map<T,T2>::iterator it = m.begin(); 
	it != m.end(); ++it) {
	v.push_back(it->first);
    }
    return v;
}

inline int hammingDistance(string & s1,string & s2){
    if(s1.size() != s2.size()){ 	cerr<<"ERROR: the hammingDistance() function cannot be called for strings of different lengths"<<endl; 	exit(1);     }
    int substitutions=0;
    for(unsigned int i=0;i<s1.size();i++){    
	if(s1[i] != s2[i])
	    substitutions++;
    }
    return substitutions;
}

inline vector<string> splitWS(string const & tosplit) { 
    istringstream tempIS(tosplit);
    vector<string> toReturn;

    copy(istream_iterator<string>(tempIS), 
	 istream_iterator<string>(),
         back_inserter(toReturn));
 
   return toReturn;
}



inline bool isStringNatNumber(string const & totest) { 

    for(unsigned int i=0;i<totest.size();i++){
	if(!isdigit(totest[i]))
	    return false;
    }
 
   return true;
}


inline bool isDos(string const & filetotest) { 
    igzstream myfiled (filetotest.c_str(),ios::in|ios::binary);

    if (!myfiled.good()){
	cerr << "Cannot open file  "<<filetotest<<""<<endl;
	exit(1);
    }

    unsigned int maxData=10000;
    unsigned int dataSoFar=0 ;

    char toread;
    char toread2=0;

    while(myfiled.read(&toread, sizeof (char))){
	//cout<<toread<<"\t"<<int(toread)<<endl;
	if(int(toread) ==  10 ){
	    if(int(toread2) == 13){
		return true;
	    }else{
		return false;
	    }
	}
	toread2 = toread;

	if( dataSoFar > maxData)
	    break;
    }
    myfiled.close();

    
    
    
    return false;
}


inline bool isMac(string const & filetotest) { 
    igzstream myfiled (filetotest.c_str(),ios::in|ios::binary);

    if (!myfiled.good()){
	cerr << "Cannot open file  "<<filetotest<<""<<endl;
	exit(1);
    }
	
    unsigned int maxData=10000;
    unsigned int dataSoFar=0 ;

    char toread;
    char toread2=0;

    while(myfiled.read(&toread, sizeof (char))){
	//cout<<toread<<"\t"<<int(toread)<<endl;
	if(int(toread2) == 13){ 
	    if(int(toread) !=  10 ){ 
		return true; 
	    }else{
		return false;
	    }
	}

	toread2 = toread;

	if( dataSoFar > maxData)
	    break;
    }
    myfiled.close();

    
    if(int(toread2) == 13)
	return true;
    else   
	return false;
}


// Returns log10( pow(10,x)+pow(10,y) ), but does so without causing
// overflow or loss of precision.
// useful if x and y are logs of probabilities
// BEWARE: of the log base, see oplusnatl()
// log( e^x + e^y)
// log( e^x (1 + e^y/e^x) )
// log( e^x ) + log (1 + e^y/e^x) )
//        x   + log (1 + e^y/e^x) )
//        x   + log1p( e^(y-x) ) 
template <typename T>
inline T oplus( T x, T y ){
    return x > y 
        ? x + log1p( pow( 10, y-x ) ) / log(10)
        : y + log1p( pow( 10, x-y ) ) / log(10) ;
}


/* // Returns log10( pow(10,x)-pow(10,y) ), but does so without causing */
/* // overflow or loss of precision. */
/* // useful if x and y are logs of probabilities */
/* // will cause error if x<y as probabilities cannot be negative */
/* template <typename T> */
/* inline T ominus( T x, T y ){ */
/*     if(x >= y)  */
/*         return (x + log1p( -1.0*pow( 10, y-x ) ) / log(10)); */
/*     else{ */
/* 	cerr<<"ERROR in ominus attempt at substracting "<<y<<" from "<<x<<" and return the log"<<endl; */
/* 	exit(1); */
/*     }	 */
/* } */



// Returns log10( pow(10,x)+pow(10,y) ), but does so without causing
// overflow or loss of precision.
// discards x if not initialized (x=0)
template <typename T>
inline T oplusInit(T x,T y ){

    if( x == 0 ){ //no initialized, as log = 0 should not exist
	return y;
    }

    return x > y 
        ? x + log1p( pow( 10, y-x ) ) / log(10)
        : y + log1p( pow( 10, x-y ) ) / log(10) ;
}


// Returns log10l( powl(10,x)+powl(10,y) ), but does so without causing
// overflow or loss of precision.
template <typename T>
inline T oplusl( T x, T y ){
    return x > y 
        ? x + log1pl( powl( 10, y-x ) ) / logl(10)
        : y + log1pl( powl( 10, x-y ) ) / logl(10) ;
}



/* // Returns log10( pow(10,x)-pow(10,y) ), but does so without causing */
/* // overflow or loss of precision. */
/* // useful if x and y are logs of probabilities */
/* // will cause error if x<y as probabilities cannot be negative */
/* template <typename T> */
/* inline T ominusl( T x, T y ){ */
/*     if(x >= y)  */
/*         return (x + log1pl( -1.0*powl( 10, y-x ) ) / logl(10)); */
/*     else{ */
/* 	cerr<<"ERROR in ominus attempt at substracting "<<y<<" from "<<x<<" and return the log"<<endl; */
/* 	exit(1); */
/*     }	 */
/* } */



// Returns log10l( powl(10,x)+powl(10,y) ), but does so without causing
// overflow or loss of precision.
// discards x if not initialized (x=0)
template <typename T>
inline T oplusInitl(T x,T y ){

    if( x == 0 ){ //no initialized, as log = 0 should not exist
	return y;
    }

    return x > y 
        ? x + log1pl( powl( 10, y-x ) ) / logl(10)
        : y + log1pl( powl( 10, x-y ) ) / logl(10) ;
}




// Returns log( expl(x)+expl(y) ), but does so without causing
// overflow or loss of precision.
template <typename T>
inline T oplusnatl( T x, T y ){
    return x > y 
        ? x + log1pl( expl( y-x ) ) 
        : y + log1pl( expl( x-y ) )  ;
}





/* // Returns log( expl(x)+expl(y) ), but does so without causing */
/* // overflow or loss of precision. */
/* template <typename T> */
/* inline T ominusnatl( T x, T y ){ */
/*     if(x > y) */
/*         return (x + log1pl( -1.0*expl( y-x ) )); */
/*     else{ */
/* 	cerr<<"ERROR in ominus attempt at substracting "<<y<<" from "<<x<<" and return the log"<<endl; */
/* 	exit(1); */
/*     } */
/* } */





// Returns log( expl(x)+expl(y) ), but does so without causing
// overflow or loss of precision.
// discards x if not initialized (x=0)
template <typename T>
inline T oplusInitnatl(T x,T y ){

    if( x == 0 ){ //no initialized, as log = 0 should not exist
	return y;
    }

    return x > y 
        ? x + log1pl( expl( y-x ) ) 
        : y + log1pl( expl( x-y ) )  ;
}




template <typename T>
inline pair<T,T> firstAndSecondHighestVector(const vector<T> & toSearch){
    if(toSearch.size() < 2 ){
	cerr<<"Cannot call firstAndSecondVector() on an array of length "<<toSearch.size()<<endl;
	exit(1);
    }

    T fstmax;
    T sncmax;

    if(toSearch[0] > toSearch[1]){
	fstmax  = toSearch[0];
	sncmax  = toSearch[1];
    }else{
	fstmax  = toSearch[1];
	sncmax  = toSearch[0];
    }
    
    for(unsigned int i=2;i<toSearch.size();i++){
	if(toSearch[i] > fstmax ){
	    sncmax = fstmax;
	    fstmax = toSearch[i];	    
	}else{
	    if(toSearch[i] >  sncmax ){
		sncmax = toSearch[i];
	    }
	}
    }

    pair<T,T> toreturn;
    toreturn.first  = fstmax;
    toreturn.second = sncmax;
    return toreturn;
}


template <typename T>
inline pair<T,T> firstAndSecondHighestArray(const T toSearch[] ,const int size){
    if(size < 2 ){
	cerr<<"Cannot call firstAndSecondArray() on an array of length "<<size<<endl;
	exit(1);
    }


    
    T fstmax;
    T sncmax;
    if(toSearch[0] > toSearch[1]){
	fstmax  = toSearch[0];
	sncmax  = toSearch[1];
    }else{
	fstmax  = toSearch[1];
	sncmax  = toSearch[0];
    }

    for(int i=2;i<size;i++){
	if(toSearch[i] > fstmax ){
	    sncmax = fstmax;
	    fstmax = toSearch[i];	    
	}else{
	    if(toSearch[i] >  sncmax ){
		sncmax = toSearch[i];
	    }
	}
    }

    pair<T,T> toreturn;
    toreturn.first  = fstmax;
    toreturn.second = sncmax;
    return toreturn;
}



inline int returnOpenFileDescriptors(){
    struct stat   statFD;

    string        toReturn="";

    int fileDMAX = getdtablesize();
     
    int currentFD = 0;
    for(int i=0;i<=fileDMAX; i++ ) {
	int statusFS=fstat(i, &statFD);

	if(statusFS == -1)
	    break;
	
	//if(errno != EBADF) 
	currentFD++;	
    }
    return currentFD;
}

inline int returnOpenFileDescriptorsMax(){

    int fileDMAX = getdtablesize();
    return fileDMAX;

}


inline string returnFileDescriptorStats(){
    struct stat   statFD;
    struct rlimit rlimitFD;
    string        toReturn="";

    int fileDMAX = getdtablesize();
     
    int currentFD = 0;
    for(int i=0;i<=fileDMAX; i++ ) {
	int statusFS=fstat(i, &statFD);

	if(statusFS == -1)
	    break;
	
	//if(errno != EBADF) 
	currentFD++;	
    }
     
    toReturn+=
	"fds currently open      :\t"+stringify(currentFD)+"\t"+
	"fds max. lim            :\t"+stringify(fileDMAX)+"\t";

    getrlimit(RLIMIT_NOFILE, &rlimitFD);

    toReturn+=
	"resource currrent limit :\t"+stringify(rlimitFD.rlim_cur) +"\t"+
	"resource max limit fds  :\t"+stringify(rlimitFD.rlim_max);
    
    return toReturn;
}

inline pair<double,double> computeJackknifeConfIntervals(double S,const vector<double> & allVals){
     pair<double,double> toreturn ;

     vector<double> Spseudo;
     double N= double(allVals.size());


     for(unsigned int i=0;i<allVals.size();i++){
	 Spseudo.push_back( S + (N-1.0)*(S-allVals[i]));
	 /* cout<<"ps\t"<<i<<"\t"<<Spseudo[i]<<"\t"<<allVals[i]<<"\tS="<<S<<"\t"<<(N-1.0)*( S-allVals[i])<<"\t"<<N<<endl; */
    }

    double pS=0;

    for(unsigned int i=0;i<allVals.size();i++){
        pS+=Spseudo[i];
    }
    pS=pS/N;
    /* cout<<"pS "<<pS<<endl; */
    double var=0;
    for(unsigned int i=0;i<allVals.size();i++){
        var+= ( pow( (Spseudo[i] - pS),2.0) );
    }
    var =var
        /
        (N-1.0) ;
    /* cout<<"var "<<var<<endl; */

    double sqrtVarDivN= sqrt(var/N);
    /* cout<<"sq var "<<sqrtVar<<endl; */
    /* cout<<"err "<<1.96*sqrtVar<<endl; */

    toreturn.first  = pS-1.96*sqrtVarDivN;
    toreturn.second = pS+1.96*sqrtVarDivN;
    /* pair<double,double> tempV=computeMeanSTDDEV(allVals); */
    /* //double sqrtVar= sqrt(var); */

    /* toreturn.first  = S-1.96*tempV.second; */
    /* toreturn.second = S+1.96*tempV.second; */

    return toreturn;
}



//taken from http://stackoverflow.com/questions/2844817/how-do-i-check-if-a-c-string-is-an-int
inline bool isInteger(const std::string & s){
   if(s.empty() || ((!isdigit(s[0])) && (s[0] != '-') && (s[0] != '+')))
       return false ;
   char * p ;
   strtol(s.c_str(), &p, 10) ;

   return (*p == 0) ;
}


//Procedure to compare chromosome names returns 
// -1 if chr1<chr2
//  0 if both chromosomes are equal
//  1 if chr1>chr2
inline int compare2Chrs(const string & chr1,const string & chr2){

    /* cout<<chr1<<"\t"<<chr2<<endl; */

    if(chr1==chr2)
	return 0;
    
    //try to convert chr1 to int

    istringstream i1(chr1);
    int           chr1int;
    bool is1anint = isInteger(chr1);
    if(is1anint)
	(i1 >> chr1int);//if worked

    istringstream i2(chr2);
    int           chr2int;
    bool is2anint = isInteger(chr2);
    if(is2anint)
	(i2 >> chr2int);//true if worked, false otherwise

    if( is1anint &&   is2anint){
	if(chr1int<chr2int)
	    return -1;
	if(chr1int>chr2int)
	    return  1;
    }

    //we imagine with 1000 Genomes, the 1..22 come before the rest MT etc
    if(!is1anint &&  is2anint){
	return   1;
    }

    //we imagine with 1000 Genomes, the 1..22 come before the rest
    if( is1anint &&  !is2anint){
	return  -1;
    }

    //we imagine with 1000 Genomes, the 1..22 come before the rest
    if( !is1anint && !is2anint){
	unsigned int commonPrefixLength=0;
	unsigned int minSize = min( chr1.size() ,chr2.size() );

	for(commonPrefixLength=0;commonPrefixLength<minSize;commonPrefixLength++){
	    if( chr1[commonPrefixLength] != chr2[commonPrefixLength] ){
		break;
	    }else{
		//just capture the common letters
		if( isdigit(chr1[commonPrefixLength] ) )
		    break;
	    }

	}
	
	if(commonPrefixLength == 0) { //no common prefix
	    if(chr1<chr2)
		return -1;
	    if(chr1>chr2)
		return  1;
	}else{
	    string chr1_ = chr1.substr(commonPrefixLength);
	    string chr2_ = chr2.substr(commonPrefixLength);
	    //cout<<chr1_<<"#"<<chr2_<<endl;

	    istringstream i1_(chr1_);
	    is1anint = isInteger(chr1_);
	    if(is1anint)
		(i1_ >> chr1int);//true if worked, false otherwise

	    istringstream i2_(chr2_);
	    is2anint = isInteger(chr2_);
	    if(is2anint)
		(i2_ >> chr2int);//true if worked, false otherwise


	    if( is1anint &&   is2anint){
		//cout<<chr1_<<"#"<<chr2_<<endl;
		/* if(chr1_<chr2_) */
		/*     return -1; */
		/* if(chr1_>chr2_) */
		/*     return  1; */
		unsigned int commonPrefixLength2=0;
		for(commonPrefixLength2=0;commonPrefixLength2<minSize;commonPrefixLength2++){
		    if( chr1_[commonPrefixLength2] != chr2_[commonPrefixLength2] ){
			break;
		    }
		    		    
		}
		//cout<<"commonPrefixLength2\t"<<commonPrefixLength2<<endl;
		

		if(commonPrefixLength2 == 0){
		    
		    if(chr1_<chr2_)
			return -1;
		    if(chr1_>chr2_)
			return  1;

		}else{
		  
		    if(chr1_.size() > chr2_.size() ){
			return -1;
		    }else{
			if(chr1_.size() < chr2_.size() ){
			    return 1;
			}else{			    
			    /* cerr<<"Invalid state#2 in compare2Chrs()"<<endl; */
			    /* exit(1); */
			    if(chr1_<chr2_)
				return -1;
			    if(chr1_>chr2_)
				return  1;
			}
		    }
		}

	    }

	    //we imagine with 1000 Genomes, the 1..22 come before the rest MT etc
	    if(!is1anint &&  is2anint){
		return   1;
	    }

	    //we imagine with 1000 Genomes, the 1..22 come before the rest
	    if( is1anint &&  !is2anint){
		return  -1;
	    }

	    //we imagine with 1000 Genomes, the 1..22 come before the rest
	    if( !is1anint && !is2anint){
		//cout<<chr1_<<"#"<<chr2_<<endl;

		if(chr1_.size() < chr2_.size() ){
		    return -1;
		}else{
		    if(chr1_.size() > chr2_.size() ){
			return 1;
		    }else{	
			if(chr1_<chr2_)
			    return -1;
			if(chr1_>chr2_)
			    return  1;
		    }
		}
	    }



	    /* cerr<<"Invalid state#1 in compare2Chrs()"<<endl; */
	    /* exit(1); */

	    //return compare2Chrs(chr1_,chr2_);
	}
    }



    cerr<<"Invalid state#1 in compare2Chrs()"<<endl;
    exit(1);

}



inline bool cmp2Chrs(const string & a,const string & b){
    int res = compare2Chrs(a,b);
    return (res == -1 );
}


inline char dinucleotide2uipac(const char & b1_,const char & b2_){
    char b1;
    char b2;

    if(b1_ < b2_){
	b1=upper(b1_);
	b2=upper(b2_);
    }else{
	b1=upper(b2_);
	b2=upper(b1_);
    }

    if(b1 == 'A'){

	if(b2 == 'A')
	    return 'A';
	if(b2 == 'C')
	    return 'M';
	if(b2 == 'G')
	    return 'R';
	if(b2 == 'T')
	    return 'W';       

    }

    if(b1 == 'C'){

	if(b2 == 'C')
	    return 'C';
	if(b2 == 'G')
	    return 'S';
	if(b2 == 'T')
	    return 'Y';       
    }

    if(b1 == 'G'){

	if(b2 == 'G')
	    return 'G';       
	if(b2 == 'T')
	    return 'K';       
    }

    if(b1 == 'T'){
	if(b2 == 'T')
	    return 'T';       
		    
    }

    cerr<<"libgab.h dinucleotide2uipac() Invalid bases "<<b1_<<" "<<b2_<<endl;
    exit(1);

}

//taken from http://www.masaers.com/2013/10/08/Implementing-Poisson-pmf.html
inline double poisson_pmf(const double k, const double lambda) {
    double lgamk = lgamma(k + 1.0);
    double lgam  = log(lambda);
    
    return exp(k*lgam-lgamk-lambda);
}

inline long double poisson_pmfl(const long double k, const long double lambda) {
    double lgamkl = lgammal(k + 1.0);
    double lgaml  = logl(lambda);
    
    return expl(k*lgaml-lgamkl-lambda);
}


// had problems
/* template <typename T> */
/* string thousandSeparator(const T i){ */
/*     stringstream s; */
/*     s.imbue(std::locale("en_US.UTF-8")); */
/*     s << i; */
/*     return s.str(); */
/* } */


template <typename T>
string thousandSeparator(const T i){
    
    if( !isInt(stringify(i)) ){
	cerr<<"ERROR: Cannot add thousandSeparator to non-integer "<<i<<endl;
	exit(1);
    }

    stringstream s;
    string       s_;
    string       sToReturn="";
    s << i; 
    s_ = s.str();

    size_t l=0;
    for(int i=int(s_.size()-1);i>=0;i--){

	sToReturn=s_[i]+sToReturn;
	if( (l%3)==2 && l!=0 && i!=0){
	    sToReturn=","+sToReturn;
	}

	l++;
    }
    
    return sToReturn;
}

//taken from https://stackoverflow.com/questions/14539867/how-to-display-a-progress-indicator-in-pure-c-c-cout-printf
inline void printprogressBarCerr(float progress){
    int barWidth = 70;

    cerr << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
	if (i < pos)
	    cerr << "=";
	else
	    if (i == pos)
		cerr << ">";
	    else
		cerr << " ";
	
    }
    cerr << "] " << int(progress * 100.0) << " %\r";
    cerr.flush();
}


//removes digits and dots at the end
inline string removeDigitsDots(const string s){
    string toreturn="";
    bool foundNondgdot=false;
    for(int i=int(s.size()-1);i>-1;i--){
	if( isdigit(s[i]) || (s[i] == '.')){
	    if(foundNondgdot) 	    toreturn = s[i]+toreturn;//we add chars once we found a non digit and dot
	}else{
	    foundNondgdot=true;
	    toreturn = s[i]+toreturn;
	}
    }

    return toreturn;
}

//ghdir contains the .git directory
inline string determineGHTag(const string  programName,const string  suffixToAdd){

    string directoryProgram;
    string commandPath=string(programName);
    size_t posSlash=commandPath.find_last_of("/");
    if(posSlash == string::npos){
        directoryProgram="";
    }else{
        directoryProgram=commandPath.substr(0,posSlash);
    }
//string gitFileLog=directoryProgram+"/"+suffixToAdd+"/.git/logs/HEAD";

    ifstream myFileLog;
    string line;
    string gitFileLog = directoryProgram+"/"+suffixToAdd+"/.git/ORIG_HEAD";
    string gitCMT;
    myFileLog.open(gitFileLog.c_str(), ios::in);
    
    if (myFileLog.is_open()){
	while ( getline (myFileLog,line)){
	    gitCMT=line;
	}
	myFileLog.close();	
    }else{
	cerr << "Unable to open github file "<<gitFileLog<<endl;
	return "NA";
    }

    ifstream myFileTag;
    string gitVersion="NA";
    string gitFileTag = directoryProgram+"/"+suffixToAdd+"/.git/packed-refs";
    myFileTag.open(gitFileTag.c_str(), ios::in);
    
    if (myFileTag.is_open()){
	while ( getline (myFileTag,line)){
	    vector<string> vs=allTokens(line,' ');
	    if(vs[0] == gitCMT){
		gitVersion=vs[1];
	    }	    
	}
	myFileTag.close();	
    }else{
	cerr << "Unable to open github file "<<gitFileTag<<endl;
	return "NA";
    }

    string rt="refs/tags/";
    if(strBeginsWith(gitVersion,rt)){
	gitVersion=gitVersion.substr(rt.length(),gitVersion.length());
    }
    return gitVersion;
}


#endif


