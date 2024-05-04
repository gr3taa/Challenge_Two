#ifndef _CHRONO_H_INCLUDE_
#define _CHRONO_H_INCLUDE_
#include <chrono>
#include <iostream>

namespace Timings
{
class Chrono
{
public:
  using MyClock = std::chrono::high_resolution_clock;
  using MyTimePoint = std::chrono::time_point<MyClock>;

private:
  MyTimePoint startTime{MyClock::now()};
  MyTimePoint stopTime{MyClock::now()};

public:
  //! Outputs time from last start and stop
  friend std::ostream &operator<<(std::ostream &, Chrono const &);
  //! Starts/reset  counting time
  inline void start();
  //! Stops counting time
  inline void stop();
  //! Outputs wall time between last start and stop (in microseconds)
  inline double wallTime() const;
  //! Outputs wall time between last start and now! (in microseconds)
  inline double wallTimeNow() const;
};
void
Chrono::start()
{
  startTime = stopTime = MyClock::now();
}

void
Chrono::stop()
{
  stopTime = MyClock::now();
}

double
Chrono::wallTime() const
{
  using namespace std::chrono;
  auto time_span = duration_cast<nanoseconds>(stopTime - startTime);
  return time_span.count() / 1000.;
}

double
Chrono::wallTimeNow() const
{
  using namespace std::chrono;
  auto time_span = duration_cast<nanoseconds>(MyClock::now() - startTime);
  return time_span.count() / 1000.;
}

inline std::ostream &
operator<<(std::ostream &out, Chrono const &c)
{
  auto oldf = out.flags();
  out << "Elapsed Time= " << std::scientific << c.wallTime() << " microsec"
      << std::endl;
  out.flags(oldf);
  return out;
}
} // namespace Timings

#endif
