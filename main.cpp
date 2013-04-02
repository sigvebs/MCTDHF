#include <src/mctdhfapplication.h>

// Library includes
#include <cstdlib>
#include <ctime>
#include <armadillo>

using namespace arma;

using namespace std;

int main(int argc, char** argv)
{
    wall_clock timer;

    MctdhfApplication *app;
    if (argc == 2)
        app = new MctdhfApplication(argv[1]);
    else
        app = new MctdhfApplication("../config.cfg");

    timer.tic();
    app->run();

    cout << "Run complete" << endl;
    cout << "Time elapsed " << timer.toc() << "s" << endl;

    delete app;
    return 0;
}
