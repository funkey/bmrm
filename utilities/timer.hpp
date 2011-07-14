/* Copyright (c) 2009, NICTA
 * All rights reserved.
 *
 * The contents of this file are subject to the Mozilla Public License Version
 * 1.1 (the "License"); you may not use this file except in compliance with
 * the License. You may obtain a copy of the License at
 * http://www.mozilla.org/MPL/
 *
 * Software distributed under the License is distributed on an "AS IS" basis,
 * WITHOUT WARRANTY OF ANY KIND, either express or implied. See the License
 * for the specific language governing rights and limitations under the
 * License.
 *
 * Authors      : Choon Hui Teo (ChoonHui.Teo@anu.edu.au)
 * Created      : 08/01/2007
 * Last Updated : 24/07/2007
 */


#ifndef _TIMER_HPP_
#define _TIMER_HPP_

#include <time.h>
#include <sys/time.h>
#include <sys/param.h>
#include <sys/times.h>
#include <sys/types.h>


#if defined(CLK_TCK)
#define TIMES_TICKS_PER_SEC double(CLK_TCK)
#elif defined(_SC_CLK_TCK)
#define TIMES_TICKS_PER_SEC double(sysconf(_SC_CLK_TCK))
#elif defined(HZ)
#define TIMES_TICKS_PER_SEC double(HZ)
#endif

/** Class for recording CPU and wall-clock time (in seconds) of program segments.
 */
class CTimer {
   private:
      timeval wallclock;
      struct tms start;
      struct tms end;

   protected:
      
      long   n;                 // number of records

      double cpu_max_time;      // longest recorded cpu time interval
      double cpu_min_time;      // shortest recorded cpu time interval
      double cpu_avg_time;      // average cpu time interval recorded
      double cpu_total_time;    // total cpu time interval recorded
      double cpu_cur_time;      // starting cpu time

      double wallclock_max_time;      // longest recorded wallclock time interval
      double wallclock_min_time;      // shortest recorded wallclock time interval
      double wallclock_avg_time;      // average wallclock time interval recorded
      double wallclock_total_time;    // total wallclock time interval recorded
      double wallclock_cur_time;      // starting wallclock time

      
   public:
      
      CTimer();             // constructor
      virtual ~CTimer(){}   // destructor
      
      void   Start();       // start recording time
      void   Stop();        // stop recording time
      void   Reset();       // reset internal time statistics

      double CPUMax();         // return cpu_max_time
      double CPUMin();         // return cpu_min_time
      double CPUAvg();         // compute and return cpu_avg_time
      double CPUTotal();       // return cpu_total_time
      double CurrentCPUTotal();// return current cpu_total_time

      double WallclockMax();         // return wall_clock_max_time
      double WallclockMin();         // return wall_clock_min_time
      double WallclockAvg();         // compute and return wallclock_avg_time
      double WallclockTotal();       // return wallclock_total_time
      double CurrentWallclockTotal();// return current wallclock_total_time
};

#endif
