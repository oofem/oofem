/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2025   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#include "progressbar.h"

#include <iostream>
#include <string>

#ifdef _WIN32
#include <windows.h>
#include <io.h>
#else
#include <sys/ioctl.h>
#include <unistd.h>
#include <signal.h>
#include <termios.h>
// struct termios termios_attrs;
#endif

namespace oofem {

ProgressBar oofem_ProgressBar;

#ifndef _WIN32
    ProgressBar* ProgressBar::active_instance = nullptr;
#endif

ProgressBar::ProgressBar() {}

void ProgressBar::initialize() {
    if (this->enabled)
        return; // already enabled

    this->enabled = true;

    if (is_stdout_redirected()) {
        // Disable pinned bar entirely
        enabled = false;
        return;
    }
 
#ifdef _WIN32
    vt_supported = enable_vt();
    if (vt_supported) {
        setup_scroll_region_win();
        hide_cursor();
    }
#else
    active_instance = this;
    // tcgetattr(STDIN_FILENO,&termios_attrs);
    signal(SIGWINCH, handle_resize);
    signal(SIGABRT,handle_crash);
    signal(SIGQUIT,handle_crash);
    signal(SIGTERM,handle_crash);
    signal(SIGINT,handle_crash);
    setup_scroll_region_unix();
    hide_cursor();
#endif
}

ProgressBar::~ProgressBar() {
    finish(); // ensure terminal state is restored exactly once
#ifndef _WIN32
    active_instance = nullptr;
    signal(SIGWINCH, SIG_DFL);
#endif
}

bool 
ProgressBar::is_stdout_redirected() {
#ifdef _WIN32
        return !_isatty(_fileno(stdout));
#else
        return !isatty(STDOUT_FILENO);
#endif
}

   
void 
ProgressBar::update(double value, const std::string& msg) {
    if (finished)
        return;

    if (!enabled)
        return; // do nothing


    if (value < 0) value = 0;
    if (value > 1) value = 1;

#ifdef _WIN32
    if (!vt_supported) {
        // Fallback: simple inline progress info, no terminal manipulation
        std::cout << msg << " " << int(value * 100) << "%\n";
        return;
    }
    draw_bar_win(value, msg);
#else
    draw_bar_unix(value, msg);
#endif
}

void ProgressBar::finish() {
    if (finished)
        return;
    if (!enabled) 
        return; // nothing to restore

#ifdef _WIN32
    if (vt_supported) {
        // restore scroll region
        std::cout << "\033[r";
        // clear last line
        std::cout << "\033[999;1H\033[2K";
        // show cursor
        show_cursor();
        std::cout.flush();
    }
#else
    // this would be the proper way:
    // tcsetattr(STDIN_FILENO,TCSANOW,&termios_attrs);

    // restore scroll region
    std::cout << "\033[r";
    // clear last line
    std::cout << "\033[999;1H\033[2K";
    // show cursor
    show_cursor();
    std::cout.flush();
#endif
    finished = true;
}
#ifdef _WIN32

    bool ProgressBar::enable_vt() {
        HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
        DWORD mode = 0;
        if (!GetConsoleMode(h, &mode))
            return false;

        DWORD newMode = mode | ENABLE_VIRTUAL_TERMINAL_PROCESSING;
        if (!SetConsoleMode(h, newMode))
            return false;

        return true;
    }

    void ProgressBar::get_size_win(int& rows, int& cols) {
        CONSOLE_SCREEN_BUFFER_INFO info;
        GetConsoleScreenBufferInfo(GetStdHandle(STD_OUTPUT_HANDLE), &info);
        rows = info.srWindow.Bottom - info.srWindow.Top + 1;
        cols = info.srWindow.Right - info.srWindow.Left + 1;
    }

    void ProgressBar::setup_scroll_region_win() {
        int rows, cols;
        get_size_win(rows, cols);

        if (rows <= 1)
            return;

        std::cout << "\033[s";              // save cursor
        std::cout << "\033[1;1H";           // move into region
        std::cout << "\033[1;" << (rows - 1) << "r";  // set region
        std::cout << "\033[u";              // restore cursor

        // If cursor was on last line, move it up
        std::cout << "\033[999D";           // col 1
        std::cout << "\033[1A";             // up one

        std::cout.flush();
    }


    void ProgressBar::draw_bar_win(double value, const std::string& msg) {
        int rows, cols;
        get_size_win(rows, cols);
        if (rows <= 0 || cols <= 0) return;

        int barWidth = cols - (int)msg.size() - 8;
        if (barWidth < 10) barWidth = 10;
        int filled = (int)(value * barWidth);

        std::cout << "\033[s";                     // save cursor
        std::cout << "\033[" << rows << ";1H";     // move to last line
        std::cout << "\033[2K";                    // clear line

        std::cout << msg << " [";
        for (int i = 0; i < barWidth; ++i)
            std::cout << (i < filled ? '#' : '-');
        std::cout << "]" << " " << int(value * 100) << "%";

        std::cout << "\033[u";                     // restore cursor
        std::cout.flush();
    }

    void ProgressBar::hide_cursor() {
        std::cout << "\033[?25l";
        std::cout.flush();
    }

    void ProgressBar::show_cursor() {
        std::cout << "\033[?25h";
        std::cout.flush();
    }
#else
    // ---------------- UNIX BACKEND ----------------

    void ProgressBar::handle_resize(int) {
        if (active_instance)
            active_instance->need_resize = true;
    }
    void ProgressBar::handle_crash(int sig) {
        if(sig!=SIGQUIT && sig!=SIGTERM && sig!=SIGINT && sig!=SIGABRT) return;
        // reset handler in case finish() crashes again for some reason
        signal(sig,SIG_DFL);
        active_instance->finish();
        // resend the signal
        raise(sig);
    }

    void ProgressBar::get_size_unix(int& rows, int& cols) {
        struct winsize ws{};
        if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &ws) == 0) {
            rows = ws.ws_row;
            cols = ws.ws_col;
        } else {
            rows = 24;
            cols = 80;
        }
    }

    void ProgressBar::setup_scroll_region_unix() {
        int rows, cols;
        get_size_unix(rows, cols);

        if (rows <= 1)
            return;

        // Save cursor
        std::cout << "\033[s";

        // Move cursor into the future scroll region
        std::cout << "\033[1;1H";

        // Set scroll region (1 .. rows-1)
        std::cout << "\033[1;" << (rows - 1) << "r";

        // Restore cursor
        std::cout << "\033[u";

        // If cursor was originally on the last line, move it up
        std::cout << "\033[999D";          // move to column 1
        std::cout << "\033[1A";            // move up one line

        std::cout.flush();
    }


    void ProgressBar::draw_bar_unix(double value, const std::string& msg) {
        if (need_resize) {
            setup_scroll_region_unix();
            need_resize = false;
        }

        const char* pctColor = 
            value < 0.33 ? "\033[31m" : // red 
            value < 0.66 ? "\033[33m" : // yellow 
            "\033[32m"; // green


        int rows, cols;
        get_size_unix(rows, cols);
        if (rows <= 0 || cols <= 0) return;

        int barWidth = cols - (int)msg.size() - 8;
        if (barWidth < 10) barWidth = 10;
        int filled = (int)(value * barWidth);

        std::cout << "\033[s";                     // save cursor
        std::cout << "\033[" << rows << ";1H";     // last line
        std::cout << "\033[2K";                    // clear

        std::cout << msg << " [";
        for (int i = 0; i < barWidth; ++i)
            std::cout << (i < filled ? "█" : "░");
        std::cout << "]" << " " << pctColor << int(value * 100) << "%\033[0m";

        std::cout << "\033[u";                     // restore cursor
        std::cout.flush();
    }

    void ProgressBar::hide_cursor() {
        std::cout << "\033[?25l";
        std::cout.flush();
    }

    void ProgressBar::show_cursor() {
        std::cout << "\033[?25h";
        std::cout.flush();
    }
#endif
};



