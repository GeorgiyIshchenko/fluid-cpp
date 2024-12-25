#include <bits/stdc++.h>
#include <vector>

using namespace std;

class FieldReader
{

public:
    FieldReader(filesystem::path path)
    {
        ifstream inputFile{ path };

        if (!inputFile.is_open())
        {
            cerr << "Error opening file" << endl;
            exit(1);
        }

        size_t n = 0;
        size_t m = 0;
        bool flag = true;

        string line;

        while (getline(inputFile, line))
        {
            // cout << " " << line << endl;
            if (flag){
                m = line.size();
                flag = false;
            }
            vector<char> row;
            for (auto&& c: line.substr(0, line.size())){
                row.push_back(c);
            }
            field.push_back(row);
            ++n;
        }

        N = n;
        M = m;
    }

public:
    size_t N;
    size_t M;

    vector<vector<char>> field;
};
