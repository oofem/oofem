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

#ifndef progressbar_h
#define progressbar_h

#include <iostream>
#include <string>
#include "oofemenv.h"


#ifdef _WIN32
//#include <windows.h>
#else
#include <sys/ioctl.h>
#include <unistd.h>
#include <signal.h>
#endif

namespace oofem {
/*
* ProgressBar provides an apt‑style pinned progress bar that stays fixed on the last terminal line while normal log output scrolls above it.
* It works on Unix terminals and modern Windows terminals (Windows Terminal, VSCode Terminal, ConPTY).
* When unsupported or when output is redirected, it automatically falls back to a safe no‑op mode.
*/
class OOFEM_EXPORT ProgressBar {
public:
    ProgressBar();

    ~ProgressBar(); 

    /**
    * Updates the pinned progress bar with the given progress value and status message.
    * If the terminal supports scroll regions (Unix, Windows Terminal, VSCode Terminal, ConPTY), the bar is rendered on the last visible line while normal log output continues scrolling above it.
    * If the terminal does not support VT sequences or if stdout is redirected, this method becomes a no‑op (or prints a simple inline fallback message, depending on configuration).
    * @param progress A floating‑point value in the range [0.0, 1.0] representing the completion percentage. Values outside the range are clamped.
    * @param statusMessage A short text label displayed before the progress bar. ANSI color escape sequences are allowed and will be rendered correctly on VT‑capable terminals. The message should end with \033[0m if color is used.
    */
    void update(double value, const std::string& msg);
    /// @brief  Restores the terminal state, removing the pinned progress bar.
    void initialize() ;
    void finish() ;
    bool is_stdout_redirected() ;

private:
    bool finished = false;
    bool enabled = false;

#ifdef _WIN32
    bool vt_supported = false;

    bool enable_vt();

    void get_size_win(int& rows, int& cols) ;

    void setup_scroll_region_win() ;


    void draw_bar_win(double value, const std::string& msg) ;

    void hide_cursor() ;

    void show_cursor() ;

#else
    // ---------------- UNIX BACKEND ----------------

    static ProgressBar* active_instance;
    bool need_resize = true;

    static void handle_resize(int) ;
    static void handle_crash(int) ;

    void get_size_unix(int& rows, int& cols) ;

    void setup_scroll_region_unix() ;

    void draw_bar_unix(double value, const std::string& msg) ;
    void hide_cursor() ;
    void show_cursor() ;
#endif
};

extern OOFEM_EXPORT ProgressBar oofem_ProgressBar;

} // namespace oofem
#endif // progressbar_h

