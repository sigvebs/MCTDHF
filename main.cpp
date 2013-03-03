
#include <src/mctdhfapplication.h>

#include <cstdlib>
#include <ctime>

using namespace std;

int main(int argc, char** argv)
{
    clock_t start = clock();

    MctdhfApplication *app;
    if (argc == 2)
        app = new MctdhfApplication(&argc, &argv, argv[1]);
    else
        app = new MctdhfApplication(&argc, &argv, "../config.cfg");

    app->run();

    cout << "Run complete" << endl;
    clock_t ends = clock();
    cout << "Time elapsed " << (double) (ends - start) / CLOCKS_PER_SEC << "s" << endl;

    delete app;
    return 0;
}
