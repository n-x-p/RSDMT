///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
///////////////////// Random Strong DMT ///////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
//// Author: Nikolai Peralta //////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

// Compile as C++14
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <assert.h>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/io.hpp>
typedef boost::numeric::ublas::compressed_matrix<bool> sparse;

// Makes printing vectors easy
template < class T >
std::ostream& operator << (std::ostream& os, const std::vector<T>& v){
    os << "<";
    for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
    {
        if(ii != v.begin()) {
            os << ", ";
        }
        os << *ii;
    }
    os << ">";
    return os;
}

// Reads input file to string
std::string fileRead(std::string filename) {
    // Reads input file to string to later be processed
    std::string facets = "";
    std::string line;
    std::string result;
    std::ifstream myfile (filename);
    if (myfile.is_open())
    {
        while ( getline (myfile,line) )
        {
            facets += line + "\n";
        }
        myfile.close();
    }
    for (int i = facets.size() - 1; i >= 0; i--) {
        if (facets[i] == '=') {
            result = facets.substr(i, facets.size() - 1);
            break;
        }
    }
    return result;
}

// Parses string from file into vector of vectors
std::vector<std::vector<int>> facetTable(std::string s) {
    // Parses string input into the facet table using the comma as a delimiter
    std::vector<std::vector<int>> table;
    s += ",";
    std::string delimiter = ",";

    size_t pos = 0;
    std::string token;
    while ((pos = s.find(delimiter)) != std::string::npos) {
        token = s.substr(0, pos);
        std::string temp = "";
        for (char& z : token){
            if (z=='[') {
                std::vector<int> x;
                table.push_back(x);
            }
            if (isdigit(z)) {
                temp += z;
            }
        }
        token = temp;
        table[table.size()-1].push_back(atoi(token.c_str()));
        s.erase(0, pos + delimiter.length());
    }
    table.erase(table.begin());
    return table;
}

// Checks to see if 2 vectors/simplices are equal
bool vEqual(std::vector<int> a, std::vector<int> b) {
    // Method for comparing vectors/simplices
    if (a.size() != b.size()) {
        return false;
    }
    for (int i = 0; i < a.size(); i++) {
        if (a[i] != b[i]) {
            return false;
        }
    }
    return true;
}

// Creates the union of 2 vectors (with redundancies)
std::vector<int> vFuse(std::vector<int> a, std::vector<int> b) {
    std::vector<int> c = {};
    for (int i = 0; i < a.size(); i++) {
        c.push_back(a[i]);
    }
    for (int i = 0; i < b.size(); i++) {
        c.push_back(b[i]);
    }
    std::sort(c.begin(),c.end());
    for (int i = 1; i < c.size(); i++) {
        if(c[i-1] == c[i]) {
            c.erase(c.begin()+i);
            i--;
        }
    }
    return c;
}

// Used for sorting vector of vectors
bool less_vectors(const std::vector<int>& a,const std::vector<int>& b) {
    return a.size() < b.size();
}

// Generates the powerset of a vector
std::vector<std::vector<int>> powerSet(std::vector<int> simplex) {
    // Returns powerset of the simplex/facet as a vector<vector<int>>
    std::vector<std::vector<int>> table = {{}};
    for (int i = 0; i < simplex.size(); i++) {
        int item = simplex[i];
        int tSize = table.size();
        for (int j = 0; j < tSize; j++) {
            std::vector<int> set = vFuse(table[j], {item});
            table.push_back(set);
        }
    }
    std::sort(table.begin(),table.end(),less_vectors);
    return table;
}

// creates union of 2 tables
std::vector<std::vector<int>> fuse(std::vector<std::vector<int>> a, std::vector<std::vector<int>> b) {
    // Adds elements of b to a
    // Essentially a Union of the sets.
    std::vector<std::vector<int>> table;
    for (int i = 0; i < a.size(); i++) {
        table.push_back(a[i]);
    }
    for (int i = 0; i < b.size(); i++) {
        table.push_back(b[i]);
    }
    std::sort(table.begin(),table.end(),less_vectors);
    for (int i = 1; i < table.size(); i++) {
        if (vEqual(table[i-1], table[i])) {
            table.erase(table.begin()+i);
            i--;
        }
    }
    return table;
}

// Checks to see if sigma > tau
bool isFace(std::vector<int> sigma, std::vector<int> tau) {
    int j;
    for (int i=0; i<tau.size(); i++)
    {
        for (j = 0; j<sigma.size(); j++)
        {
            if(tau[i] == sigma[j]) {
                break;
            }
        }
        if (j == sigma.size()) {
            return false;
        }
    }
    return true;
}

// Gives face vector of given complex/table
std::vector<int> faceVector(std::vector<std::vector<int>> table) {
    int size = 0;
    for (int i = 0; i < table.size(); i++) {
        if (table[i].size() > size) {
            size = table[i].size();
        }
    }
    std::vector<int> fv (size,0);
    for (int i = 0; i < table.size(); i++) {
        fv[table[i].size()-1]++;
    }
    return fv;
}

// Checks for redundancies after all powersets are generted
std::vector<std::vector<int>> redundancyCheck(std::vector<std::vector<int>> table) {
    // Removes redundancies in the vector of simplices
    std::vector<std::vector<int>> output = {{}};
    bool found;
    for (int i = 0; i < table.size(); i++) {
        found = false;
        for (int j = 0; j<output.size(); j++) {
            if (vEqual(table[i], output[j])) {
                found = true;
            }
        }
        if (!found) {
            output.push_back(table[i]);
        }
    }
    return output;
}

// Creates a Hasse diagram/ sparse matrix that shows what's connected to what
sparse genHasse(std::vector<std::vector<int>> table) {
    boost::numeric::ublas::compressed_matrix<bool> m (table.size(), table.size(), 0);
    for (int i = 0; i < m.size1(); i++) {
        for (int j = 0; j < m.size2(); j++) {
            if (table[j].size() == table[i].size() + 1 && isFace(table[j],table[i])) {
                m(i,j) = true;
            }
        }
    }
    return m;
}

// returns true if a&b are a free pair
bool freePair(sparse G, int a, int b) {
    //a is a subset of b
    for (int i = 0; i < G.size2(); i++){
        if (G(b,i)) {
            return false;
        }
        if (i == b) {
            continue;
        }
        if (G(a,i)) {
            return false;
        }
    }
    return true;
}

// returns true if s is a facet
bool isFacet(sparse m, int s) {
    for (int i = 0; i < m.size2(); i++) {
        if (m(s,i)) {
            return false;
        }
    }
    for (int i = 0; i < s; i++) {
        if(m(i,s)) {
            return true;
        }
    }
    return false;
}

// returns table of indices of free pairs
std::vector<std::vector<int>> fpTable (std::vector<std::vector<int>> T, sparse m) {
    std::vector<std::vector<int>> out = {};
    for (int i = 0; i < m.size1(); i++) {
        for (int j = 0; j < m.size2(); j++) {
            if (T[i].size() != T[j].size() - 1) {
                continue;
            }
            if (freePair(m,i,j)) {
                out.push_back({i,j});
            }
        }
    }
    return out;
}

// returns vector that points to indices of facets in hasse diagram
std::vector<int> facetIndex(sparse m) {
    std::vector<int> indices = {};
    for (int i = 0; i < m.size1(); i++) {
        if (isFacet(m, i)) {
            indices.push_back(i);
        }
    }
    return indices;
}

// returns the vertex set
std::vector<std::vector<int>> vertexSet(std::vector<std::vector<int>> T) {
    std::vector<std::vector<int>> vs = {};
    for (int i = 0; i < T.size(); i++) {
        if (T[i].size() == 1) {
            vs.push_back(T[i]);
        }
    }
    return vs;
}

// returns true if x dominates y
bool xDOMy(std::vector<std::vector<int>> T, sparse m, int x, int y) {
    //checks to see if x dominates y
    std::vector<std::vector<int>> vs = vertexSet(T);
    std::vector<int> fctindx = facetIndex(m);
    std::vector<int> xArr = {};
    std::vector<int> yArr = {};
    bool t1 = false;
    bool t2 = false;
    for (int i = 0; i < m.size1(); i++) {
        if (m(x,i)) {
            t1 = true;
        }
        if (m(y,i)) {
            t2 = true;
        }
    }
    if (!(t2 && t1)){
        return false;
    }
    for (int i = 0; i <fctindx.size(); i++) {
        for (int j = 0; j < T[fctindx[i]].size(); j++) {
            if (vs[x][0] == T[fctindx[i]][j]) {
                xArr.push_back(fctindx[i]);
            }
            if (vs[y][0] == T[fctindx[i]][j]) {
                yArr.push_back(fctindx[i]);
            }
        }
    }
    return isFace(xArr,yArr)||vEqual(xArr,yArr);
}

// performs strong collapse
sparse sCollapse(std::vector<std::vector<int>> T, sparse m, int v) {
    std::vector<int> oStarIndex = {};
    sparse newAM = m;
    for (int i = 0; i < T.size(); i++) {
        if (isFace(T[i], T[v])) {
            oStarIndex.push_back(i);
        }
    }
    //std::cout<<oStarIndex<<std::endl;
    for (int i = 0; i < oStarIndex.size(); i++) {
        for (int j = 0; j < newAM.size2(); j++) {
            newAM(j,oStarIndex[i]) = false;
        }
    }
    return newAM;
}

// performs elementary collapse
sparse eCollapse(sparse m, int s, int t) {
    sparse newAM = m;
    for (int i = 0; i < newAM.size1(); i++) {
        newAM(i,s) = false;
        newAM(i,t) = false;
    }
    return newAM;
}

// removes random facet
sparse lastResort(std::vector<std::vector<int>> T, sparse m, std::vector<int> * dmv) {
    //removes random facet
    //FIXME: Random top face
    srand (time(NULL));
    sparse newAM = m;
    std::vector<int> fidx = facetIndex(newAM);
    int rnum = rand() % fidx.size();
    for (int i = 0; i < newAM.size1(); i++) {
        newAM(i, fidx[rnum]) = false;
    }
    return newAM;
}

// checks to see if the sparse matrix is empty/ aka the complex is collapsed
bool isCollapsed(sparse m) {
    for (int i = 0; i < m.size1(); i++) {
        for (int j = 0; j < m.size2(); j++) {
            if (m(i,j)) {
                return false;
            }
        }
    }
    return true;
}

// returns table of dom vertices indices
std::vector<std::vector<int>> domTable(sparse m, std::vector<std::vector<int>> T) {
    std::vector<std::vector<int>> L = {};
    std::vector<std::vector<int>> vs = {};
    for (int i = 0; i < T.size(); i++) {
        if(T[i].size() == 1) {
            vs.push_back(T[i]);
        }
    }
    for (int i = 0; i < vs.size(); i++) {
        for (int j = 0; j < vs.size(); j++) {
            if (i == j) {
                continue;
            }
            if (xDOMy(T, m, i, j)) {
                L.push_back({i,j});
            }
        }
    }
    return L;
}

// returns readable table of dom vertices indices
std::vector<std::vector<int>> domTableR(sparse m, std::vector<std::vector<int>> T) {
    std::vector<std::vector<int>> L = domTable(m, T);
    for (int i = 0; i < L.size(); i++) {
        L[i][0] = T[L[i][0]][0];
        L[i][1] = T[L[i][1]][0];
    }
    return L;
}

// returns dimension of simplex
int  dim(std::vector<std::vector<int>> T) {
    int top = T[T.size()-1].size() - 1;
}

// MAIN
int main() {
    std::string filename = "scTest"; //"Barnette_sphere";
    std::string facets = fileRead("LoT/"+filename);
    std::vector<std::vector<int>> table = facetTable(facets);
    for (int i = 0;i < table.size();i++) {
        table = fuse(table,powerSet(table[i]));
    }
    table = redundancyCheck(table);
    //printTable("Full Complex", table);
    table.erase(table.begin());
    std::cout<<std::endl<<"Current face vector for the "<<filename<<" = "
             << faceVector(table) << std::endl;
    sparse hasse = genHasse(table);
    int dimension = dim(table);
    std::vector<int> DMV(dim(table) + 1, 0);
    srand (time(NULL));
    int rnum;
    std::vector<std::vector<int>> dt;
    std::vector<std::vector<int>> ft;
    std::cout<<table<<std::endl;

    //std::cout<<hasse<<std::endl;
    std::cout<<"\ndomTable: "<<domTableR(hasse,table)<<std::endl;
    std::cout<<"\nFacet Index: "<<facetIndex(hasse)<<std::endl;
    //hasse = sCollapse(table,hasse,3);
    //std::cout<<"Strong collapse on vertex "<<table[3]<<std::endl;
    //std::cout<<hasse<<std::endl;
    //std::cout<<"\ndomTable: "<<domTable(hasse,table)<<facetIndex(hasse)<<std::endl;
    std::cout<<"\nFacets: ["<<facetIndex(hasse).size()<<"] ";
    for (int i =0; i<facetIndex(hasse).size();i++){
        std::cout<<table[facetIndex(hasse)[i]];
    }
    //Actual algorithm
    // Which is broken
    /*
    while (!isCollapsed(hasse)) {
        dt = domTable(hasse, table);
        ft = fpTable(table, hasse);
        if (dt.size() > 0) {
            rnum = rand() % dt.size();
            hasse = sCollapse(table,hasse,dt[rnum][1]);
            std::cout<<"Strong collapse on vertex "<<table[dt[rnum][1]]<<std::endl;
        } else if (ft.size() > 0) {
            rnum = rand() % ft.size();
            hasse = eCollapse(hasse, ft[rnum][0], ft[rnum][1]);
            //std::cout<<"Elementary collapse with free pair "<<table[ft[rnum][0]]<<
            //
            //         ", "<<table[ft[rnum][1]]<<std::endl;//<<ft<<std::endl;
            DMV[table[ft[rnum][1].size() - 1] = 1 + DMV[table[ft[rnum][1]].size() - 1]
        } else {
            std::vector<int> fidx = facetIndex(hasse);
            rnum = rand() % fidx.size();
            for (int i = 0; i < hasse.size1(); i++) {
                hasse(i, fidx[rnum]) = false;
            }
            DMV[table[fidx[rnum]].size() - 1] = 1 + DMV[table[fidx[rnum]].size() - 1];
            std::cout<<"Removal of facet "<<table[fidx[rnum]]<<std::endl;
        }
    }
     */
    //std::cout<<"\nDiscrete Morse vector = "<<DMV<<std::endl;
    return 0;
}
