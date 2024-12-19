#include <bits/stdc++.h>

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

        vector<vector<char>> matrix;

        string line;

        while (getline(inputFile, line))
        {
            //cout << " " << line << endl;
            vector<char> row;
            stringstream ss(line);
            char num;

            while (ss >> num)
            {
                row.push_back(num);
                if (flag)
                {
                    m++;
                }
            }
            matrix.push_back(row);
            n++;
            flag = false;
        }

        N = n;
        M = m;

        field = matrix;
    }

public:
    size_t N;
    size_t M;

    vector<vector<char>> field;
};
