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
 * Last Updated : 10/01/2009
 */


#ifndef _TIMER_CPP_
#define _TIMER_CPP_

#include "timer.hpp"
#include <stdio.h>


CTimer::CTimer()
   :n(0),
    cpu_max_time(0),
    cpu_min_time(99999999),
    cpu_avg_time(0),
    cpu_total_time(0),
    cpu_cur_time(0),
    wallclock_max_time(0),
    wallclock_min_time(99999999),
    wallclock_avg_time(0),
    wallclock_total_time(0),
    wallclock_cur_time(0)
{}


// void CTimer::Start()
// {
//     cpu_cur_time = (double)clock()/CLOCKS_PER_SEC;
//     gettimeofday(&wallclock, NULL);
//     wallclock_cur_time = double(wallclock.tv_sec) + double(wallclock.tv_usec)/1e6;
// }

// void CTimer::Stop()
// {
//     double stopTime = (double)clock()/CLOCKS_PER_SEC;
//     cpu_cur_time = stopTime - cpu_cur_time;
//     gettimeofday(&wallclock, NULL);
//     wallclock_cur_time = (double(wallclock.tv_sec) + double(wallclock.tv_usec)/1e6 - wallclock_cur_time);
    
//     n++;
    
//     cpu_total_time += cpu_cur_time;
//     if(cpu_cur_time > cpu_max_time) cpu_max_time = cpu_cur_time;
//     if(cpu_cur_time < cpu_min_time) cpu_min_time = cpu_cur_time;
    
//     wallclock_total_time += wallclock_cur_time;
//     if(wallclock_cur_time > wallclock_max_time) wallclock_max_time = wallclock_cur_time;
//     if(wallclock_cur_time < wallclock_min_time) wallclock_min_time = wallclock_cur_time;
//     }



void CTimer::Start()
{   
   times(&start); 
   cpu_cur_time = (double(start.tms_utime) + double(start.tms_stime))/TIMES_TICKS_PER_SEC;
   gettimeofday(&wallclock, NULL);
   wallclock_cur_time = double(wallclock.tv_sec) + double(wallclock.tv_usec)/1e6;
     
}


void CTimer::Stop()
{  
   times(&end);
   cpu_cur_time = (double(end.tms_utime) + double(end.tms_stime))/TIMES_TICKS_PER_SEC - cpu_cur_time;
   
   gettimeofday(&wallclock, NULL);
   wallclock_cur_time = (double(wallclock.tv_sec) + double(wallclock.tv_usec)/1e6 - wallclock_cur_time);
   
   n++;
   
   cpu_total_time += cpu_cur_time;
   if(cpu_cur_time > cpu_max_time) cpu_max_time = cpu_cur_time;
   if(cpu_cur_time < cpu_min_time) cpu_min_time = cpu_cur_time;
   
   wallclock_total_time += wallclock_cur_time;
   if(wallclock_cur_time > wallclock_max_time) wallclock_max_time = wallclock_cur_time;
   if(wallclock_cur_time < wallclock_min_time) wallclock_min_time = wallclock_cur_time;
   
}


void CTimer::Reset()
{
   n = 0;
   
   cpu_max_time   = 0;
   cpu_min_time   = 99999999;
   cpu_avg_time   = 0;
   cpu_total_time = 0;
   cpu_cur_time   = 0;
   
   wallclock_max_time   = 0;
   wallclock_min_time   = 99999999;
   wallclock_avg_time   = 0;
   wallclock_total_time = 0;
   wallclock_cur_time   = 0;
   
}


double CTimer::CPUMax()
{
   return cpu_max_time;
}

 
double CTimer::CPUMin()
{
   return cpu_min_time;
}


double CTimer::CPUAvg()
{
   return cpu_total_time/n;
}                


double CTimer::CPUTotal()
{
   return cpu_total_time;
}


double CTimer::CurrentCPUTotal()
{
   times(&end);
   double time = (double(end.tms_utime) + double(end.tms_stime))/TIMES_TICKS_PER_SEC - cpu_cur_time;
   return time;
}



double CTimer::WallclockMax()
{
   return wallclock_max_time;
}
 

double CTimer::WallclockMin()
{
   return wallclock_min_time;
}


double CTimer::WallclockAvg()
{
   return wallclock_total_time/n;
}                


double CTimer::WallclockTotal()
{
   return wallclock_total_time;
}

double CTimer::CurrentWallclockTotal()
{
   gettimeofday(&wallclock, NULL);
   double time = (double(wallclock.tv_sec) + double(wallclock.tv_usec)/1e6 - wallclock_cur_time);
   return time;
}

#endif
