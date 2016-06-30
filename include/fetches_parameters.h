#ifndef FETCHES_PARAMETERS_H_
#define FETCHES_PARAMETERS_H_

#include <string>
#include <vector>
#include <sqlite3.h>
#include "Model.h"


class FetchesParameters {
    public:
        FetchesParameters(const std::string &filename);
        ~FetchesParameters();
        std::vector<Model> fetch_models();

    private:
        Model model_from_statement(sqlite3_stmt *stmt);

        std::string filename;
        sqlite3 *ppDb;

};



#endif //  FETCHES_PARAMETERS_H_
