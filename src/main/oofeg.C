/* $Header: /home/cvs/bp/oofem/main/src/oofeg.C,v 1.11.4.1 2004/04/05 15:19:41 bp Exp $ */
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
 *               Copyright (C) 1993 - 2012   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

//  MAIN
//  Solves finite element problems.
//
//  The DEBUG option MUST be used (check in file 'debug.def').
//  See also file 'main2.c'.
#ifdef __OOFEG

#include "oofemtxtdatareader.h"
#include "engngm.h"
#include "timestep.h"
#include "freestor.h"
#include "compiler.h"
#include "error.h"
#include "oofeggraphiccontext.h"


#include "conTable.h"
#include "mathfem.h"
#include "zznodalrecoverymodel.h"
#include "nodalaveragingrecoverymodel.h"
#include "sprnodalrecoverymodel.h"

#include "util.h"
#include "node.h"
#include "intarray.h"
#include "oofemdef.h"
#include "errorestimator.h"
#include "remeshingcrit.h"
#include "drawmode.h"
#include "contextioerr.h"
#include "oofem_terminate.h"
#include "classfactory.h"
//
// for c++ compiler to be succesfull on some c files
//

//#define OOFEG_DEVEL

using namespace oofem;

#define OOFEG_BACKGROUND_COLOR "WhiteSmoke"
#define OOFEG_DEFAULTDRAW_COLOR "black"

#define OOFEG_PALETTE_BUTTON_RESOURCE "oofeg_palette"
#define OOFEG_MENU_BUTTON_RESOURCE    "oofeg_menu"

#define OOFEG_WIRED_RENDER 0
#define OOFEG_NORMAL_RENDER 1
#define OOFEG_HIDDEN_RENDER 2
#define OOFEG_SHADED_RENDER 3


static int oofeg_wired_render = OOFEG_WIRED_RENDER;
static int oofeg_normal_render = OOFEG_NORMAL_RENDER;
static int oofeg_hidden_render = OOFEG_HIDDEN_RENDER;
static int oofeg_shaded_render = OOFEG_SHADED_RENDER;

static Widget deformations_menu, activeStep_palette,
              layers_palette, plot_palette, varplot_palette, filters_palette,
              matregfilter_palette, dofman_menu,
              color_scale_setup_palette, grey_scale_setup_palette, scale_min, scale_max,
              greyscale_min, greyscale_max,
              scale_setup_ok, toggle_scale_palette,
              animate_setup_palette, start_step, end_step, animate_scale_setup_ok,
              smoother_palette, stressForce_palette, bendingMoment_palette,
              scalarplot_palette, vectorplot_palette, tensorplot_palette,
              scalarstrain_palette, plotalgo_palette, valmode_palette, scalarerror_palette,
              scalardamage_palette, scalarplasticstrain_palette, color_scale_palette,
              oofeg_file_palette, frame_palette, view_palette, render_palette, bgcolor_palette,
              layers_update_palette;
#ifdef OOFEG_DEVEL
static Widget debug_palette;
#endif

static void OOFEGReturnHitInCmd(Widget w, XEvent *event, String *params,
                                Cardinal *num_params);

static XtTranslations tt1;
static XtActionsRec oofeg_remap_return[] = {
    { ( String ) "oofegretActCmd", OOFEGReturnHitInCmd },
};

static int oofeg_box_setup = 0;
static int oofeg_axes      = 1;

EngngModel *problem;



void deleteGraphics(oofegGraphicContext &gc);
void GeoPlot(Widget wid, XtPointer cl, XtPointer cd);
//void AugmentCommandTable();
void defPlot(Widget wid, XtPointer cl, XtPointer cd);
void  defAutoScale(Widget wid, XtPointer cl, XtPointer cd);
void  setNumerOfContours(Widget wid, XtPointer cl, XtPointer cd);
void eigVecPlot(Widget wid, XtPointer cl, XtPointer cd);
void nodeAnnPlot(Widget wid, XtPointer cl, XtPointer cd);
void elementAnnPlot(Widget wid, XtPointer cl, XtPointer cd);
void nodePlot(Widget wid, XtPointer cl, XtPointer cd);
void bcPlot(Widget wid, XtPointer cl, XtPointer cd);
//void stressPlot (int);
void varPlot(Widget w, XtPointer ptr, XtPointer call_data);
void beamForcePlot(Widget w, XtPointer ptr, XtPointer call_data);
void stresscompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void princstresscompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void equivstraincompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void straincompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void princstraincompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void plaststraincompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void princplaststraincompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void damagecompPlot(Widget w, XtPointer ptr, XtPointer call_data);
//mj
void epseqcompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void kappacompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void kappa2compPlot(Widget w, XtPointer ptr, XtPointer call_data);
void dissworkcompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void stressworkcompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void freeenergycompPlot(Widget w, XtPointer ptr, XtPointer call_data);
//emj
void princdamagecompPlot(Widget w, XtPointer ptr, XtPointer call_data);
void emptyaction(Widget w, XtPointer ptr, XtPointer call_data);
void plotAlgosel(Widget w, XtPointer ptr, XtPointer call_data);
void valmodesel(Widget w, XtPointer ptr, XtPointer call_data);
void print_state(Widget w, XtPointer text_ptr, XtPointer call_data);
void do_print_state(EView *v_p, caddr_t data_p);
void OOFEGSimpleCmd(char *);
static void apply_change(Widget w, XtPointer ptr, XtPointer call_data);
static void apply_layer_settings(Widget w, XtPointer ptr, XtPointer call_data);
static void apply_layer_update(Widget w, XtPointer ptr, XtPointer call_data);
static void apply_mat_reg_filter(Widget w, XtPointer ptr, XtPointer call_data);
char *readSimpleString(char *source, char *simpleString, char **remain);
void deleteLayerGraphics(int iLayer);
void nextStep(Widget wid, XtPointer cl, XtPointer cd);
void previousStep(Widget wid, XtPointer cl, XtPointer cd);
static void set_layer_on_off(EView *v_p, caddr_t data, WCRec *p);
static void uninstall_apply_to_view(EView *v_p, caddr_t data);
static void pass_setscale_command(Widget w, XtPointer ptr, XtPointer call_data);
static void pass_setgreyscale_command(Widget w, XtPointer ptr, XtPointer call_data);
static void pass_setanimate_command(Widget w, XtPointer ptr, XtPointer call_data);
static void apply_toggleIcs(Widget w, XtPointer ptr, XtPointer call_data);
static void toggleIcs(EView *v_p, caddr_t data, WCRec *p);
static void toggle_IcsColors(Widget w, XtPointer ptr, XtPointer call_data);
static void setIcsToColor(Widget w, XtPointer ptr, XtPointer call_data);
static void setIcsToGrey(Widget w, XtPointer ptr, XtPointer call_data);
static void toggleSmoothScale(Widget w, XtPointer ptr, XtPointer call_data);
static void toggleTransparentContours(Widget w, XtPointer ptr, XtPointer call_data);
static void setAutoscale(Widget w, XtPointer ptr, XtPointer call_data);
static void set_bg_color(Widget w, XtPointer ptr, XtPointer call_data);
static int set_background(NODE _data, NODE _view);
static void toggleStatus(Widget wid, XtPointer cl, XtPointer cd);
static void view_toggleStatus_on_off(EView *v_p, caddr_t data, WCRec *p);
#ifdef OOFEG_DEVEL
void debug_run(Widget w, XtPointer ptr, XtPointer call_data);
#endif

void setSmoother(Widget w, XtPointer ptr, XtPointer call_data);
void setSmoother(SmootherType mode);
void setupSmoother(oofegGraphicContext &gc);
void setupData(oofegGraphicContext &gc);
void drawData(oofegGraphicContext &gc);
int  updateDefPlotFlag();
void showSparseMtrxStructure(Widget wid, XtPointer cl, XtPointer cd);
void errorcompPlot(Widget w, XtPointer ptr, XtPointer call_data);

static Widget oofeg_add_palette(const char *palette_label, Widget parent, Widget *palette);
static Widget oofeg_add_palette(const char *palette_label, Widget parent, Widget *palette);
static Widget oofeg_add_popdown_menu(const char *menu_label, Widget parent, Widget *palette);
static Widget oofeg_add_button(const char *name, const char *button_label, WidgetClass wclass, Widget palette,
                               XtCallbackProc action, XtPointer data);
static Widget oofeg_add_button(const char *name, const char *button_label, WidgetClass wclass, Widget palette, Arg *arg, int ac,
                               XtCallbackProc action, XtPointer data);
static Widget oofeg_add_menu_item(const char *name, const char *item_label, Widget palette,
                                  XtCallbackProc action, XtPointer data);
static Widget oofeg_add_dialog(const char *name, const char *dialog_label, const char *prompt, const char *init_value, Widget palette,
                               XtCallbackProc action, const char *data, ESIDialogValueType type, ESIVerifyValueProc proc);

void oofeg_display_message(const char *message);
void oofeg_exit(Widget w, XtPointer ptr, XtPointer call_data);
void oofeg_quit(XtPointer ptr);
void oofeg_open_frame(Widget w, XtPointer ptr, XtPointer call_data);
void oofeg_close_frame(Widget w, XtPointer ptr, XtPointer call_data);
void oofeg_destroy_frame(EView *view, caddr_t data, WCRec *p);
void oofeg_set_render(Widget w, XtPointer ptr, XtPointer call_data);
void oofeg_toggle_axes(Widget w, XtPointer ptr, XtPointer call_data);
void set_render(EView *v_p, caddr_t data, WCRec *p);
int show_axes(NODE _data, NODE _view);
int hide_axes(NODE _data, NODE _view);
void oofeg_fit_all_graphics(Widget w, XtPointer ptr, XtPointer call_data);
void updateISA(oofegGraphicContext *context);
void updateGraphics();

static const char *OOFEG_layer_names[] = {
    "GEOMETRY LAYER", "DEFORMED GEOMETRY", "NODE ANNOTATION", "ELEMENT ANNOTATION", "VARPLOT CONTOURS", "CRACK PATTERNS", "IC-BC ANNOTATIONS", "NATURAL_BC", "SPARSE PROFILE LAYER", "DEBUG LAYER"
};


static SmootherType oofeg_smoother_modes [] = {
    Smoother_NA, Smoother_ZZ, Smoother_SPR
};


static int vectorAddr[] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12
};


#define OOFEG_BG_COLORS 3

static char oofeg_bg_color_name [ OOFEG_BG_COLORS ] [ 32 ] = {
    "white",
    "darkslategrey",
    "midnightblue"
};



static EPixel oofeg_background_color;

static DrawMode oofeg_draw_modes[] = {
    sxForce,  syForce,  szForce, // 0 - 2
    syzForce, szxForce, sxyForce, // 3 - 5

    mxForce, myForce, mzForce, // 6 - 8
    myzForce, mzxForce, mxyForce, // 9 - 11

    yieldState, crackedState,  // 12 - 13
    damageLevel,               // 14
    errorIndicatorLevel,       // 15
    relativeMeshSizeDensity,   // 16
    temperatureField,          // 17
    massConcentration1Field,    // 18

    velocityField,              // 19
    pressureField,              // 20
    vofField,                   // 21
    densityField                // 22
};

std::string jobName = "";
std::string viewTitle = "";

oofegGraphicContext gc [ OOFEG_LAST_LAYER ];
EView *myview;

/* Defaul oofem loggers */
Logger oofem :: oofem_logger(Logger :: LOG_LEVEL_INFO, stdout);
Logger oofem :: oofem_errLogger(Logger :: LOG_LEVEL_WARNING, stderr);

/* Global class factory */
ClassFactory oofem :: classFactory;

int
main(int argc, char *argv[])
{
    int i;
    int inputFileFlag = 0;
    //int rank=0;
    std :: stringstream inputFileName;
    char buff [ 20 ];
    unsigned long mask = 0;
    bool parallelFlag


    // print prg header on stdout
    printf("%s", PRG_HEADER_SM);

#ifdef __PARALLEL_MODE
    parallelFlag = true; ///@todo Read this from input arguments (eventually)
    int rank = 0;
 #ifdef __USE_MPI
    MPI_Init(& argc, & argv);
    MPI_Comm_rank(MPI_COMM_WORLD, & rank);
 #endif
#endif

#ifdef __PETSC_MODULE
    PetscInitialize(& argc, & argv, PETSC_NULL, PETSC_NULL);
#endif

    //
    // check for options
    //
    if ( argc != 1 ) {
        for ( i = 1; i < argc; i++ ) {
            if ( strcmp(argv [ i ], "-f") == 0 ) {
                if ( i + 1 < argc ) {
                    inputFileName << argv [ i + 1 ];
                    inputFileFlag = 1;
                }
            } else if ( strcmp(argv [ i ], "-l") == 0 ) {
                if ( i + 1 < argc ) {
                    strcpy(buff, argv [ i + 1 ]);
                    int level = strtol(buff, ( char ** ) NULL, 10);
                    oofem_logger.setLogLevel(level);
                    oofem_errLogger.setLogLevel(level);
                }
            }
        }
    }

    if ( !inputFileFlag ) {
        std::string input;
        :: giveInputDataFileName(input);
        inputFileName << input;
    }

#ifdef __PARALLEL_MODE
    if ( parallelFlag ) {
        inputFileName << "." << rank;
    }
#endif

    OOFEMTXTDataReader dr(inputFileName.str().c_str());

    // extract job name
    std::string temp = inputFileName.str();
    // extract job name
    size_t k = temp.size();
    while ( k-- ) {
        if ( temp[k] == '/' ) {
            k++;
            break;
        }
    }
    std::string jobName = temp.substr(k);
    std::string viewTitle = "OOFEG ("+jobName+")";

    problem = InstanciateProblem(& dr, _postProcessor, 0, parallelFlag);
    dr.finish();
    problem->checkProblemConsistency();

    mask = ESI_GRAPHIC_EDITOR_MASK;
#ifndef OOFEG_DEVEL
    mask = 0;
    mask |= ESI_TRACK_AREA_MASK;
    mask |= ESI_PROMPT_AREA_MASK;
    mask |= ESI_STATUS_AREA_MASK;
    /* mask |= ESI_COMMAND_AREA_MASK; */
    mask |= ESI_CUSTOM_PALETTE_BOX_MASK;

    oofeg_box_setup = 1;
#endif



    ESIBuildInterface(mask, argc, argv);
    myview  =  ElixirNewView(viewTitle.c_str(), const_cast< char * >("OOFEG"), const_cast< char * >(OOFEG_BACKGROUND_COLOR),
                             const_cast< char * >(OOFEG_DEFAULTDRAW_COLOR), 500, 400);
    EVSetRenderMode(myview, WIRE_RENDERING);
    EMAttachView(age_model, myview);
    gc [ 0 ].init(problem); // init all gcs

    // AugmentCommandTable();

    EMSetAssocFringeTable( age_model, gc [ 0 ].getFringeTable() );
    EVSetAssocFringeTable( myview, gc [ 0 ].getFringeTable() );

    setSmoother(Smoother_ZZ);
    //gc.setSmoother (new ZZNodalRecoveryModel(problem));
    //problem -> instanciateYourself () ;
    //problem -> giveEngngModel()->forceEquationNumbering();

    // activate some useful layers
    for ( i = 0; i < OOFEG_LAST_LAYER; i++ ) {
        EVSetLayerOnOff(myview, i, 1);
    }

    if ( oofeg_box_setup ) {
        oofeg_display_message( oofem_tmpstr(OOFEG_VERSION) );
    }

    updateISA(gc);

    //ESIPopupAndRun();

    // this versiont supports keybort shortcuts like Ctrl+A, Ctrl+X etc.
    ESIPopup();
    ESIEventLoop(true, NULL);

#ifdef __PETSC_MODULE
    PetscFinalize();
#endif
#ifdef __PARALLEL_MODE
    MPI_Finalize();
#endif

    return 0;
}

/*
 * static top_command   table[] = {
 * };
 *
 * void AugmentCommandTable(void){
 * TypeInAugmentCommandTable(table, XtNumber(table));
 * }
 */

void ESICustomize(Widget parent_pane)
{
    char tname [ 32 ], ltname [ 32 ];
    int i, nmat;

    if ( oofeg_box_setup ) {
        oofeg_add_palette("< File >", parent_pane, & oofeg_file_palette);
        oofeg_add_button("OOFEG_EXIT", "Exit", commandWidgetClass, oofeg_file_palette, oofeg_exit, ( XtPointer ) parent_pane);
    }

    oofeg_add_palette("< View >", parent_pane, & view_palette);
    {
        oofeg_add_palette("< Frame >", view_palette, & frame_palette);
        oofeg_add_button("OPEN_FRAME", "Open frame", commandWidgetClass, frame_palette, oofeg_open_frame, NULL);
        oofeg_add_button("CLOSE_FRAME", "Close frame", commandWidgetClass, frame_palette, oofeg_close_frame, NULL);

        oofeg_add_button("FIT_ALL", "Fit all", commandWidgetClass, frame_palette, oofeg_fit_all_graphics, NULL);
        oofeg_add_palette("< Render >", view_palette, & render_palette);
        oofeg_add_button("RENDER_BUTTON", "Wired", commandWidgetClass, render_palette,
                         oofeg_set_render, ( XtPointer ) & oofeg_wired_render);
        oofeg_add_button("RENDER_BUTTON", "Normal", commandWidgetClass, render_palette,
                         oofeg_set_render, ( XtPointer ) & oofeg_normal_render);
        oofeg_add_button("RENDER_BUTTON", "Hidden", commandWidgetClass, render_palette,
                         oofeg_set_render, ( XtPointer ) & oofeg_hidden_render);
        oofeg_add_button("RENDER_BUTTON", "Shaded", commandWidgetClass, render_palette,
                         oofeg_set_render, ( XtPointer ) & oofeg_shaded_render);

        oofeg_add_button("AXES_BUTTON", "Axes", commandWidgetClass, view_palette, oofeg_toggle_axes, NULL);
        oofeg_add_popdown_menu("< Bg Color >", view_palette, & bgcolor_palette);
        {
            for ( i = 0; i < OOFEG_BG_COLORS; i++ ) {
                oofeg_add_menu_item(NULL, oofeg_bg_color_name [ i ], bgcolor_palette, set_bg_color, ( XtPointer ) oofeg_bg_color_name [ i ]);
            }
        }

        oofeg_add_palette("< Color scale >", view_palette, & color_scale_palette);
        {
            oofeg_add_button("TOGGLE_SCALE", "Toggle scale",
                             commandWidgetClass, color_scale_palette, apply_toggleIcs, ( XtPointer ) toggle_scale_palette);

            oofeg_add_button("TOGGLE_COLOR_SCALE", "Revert scale colors",
                             commandWidgetClass, color_scale_palette, toggle_IcsColors, ( XtPointer ) NULL);

            oofeg_add_button("TOGGLE_SCALE_COLOR", "Color scale",
                             commandWidgetClass, color_scale_palette, setIcsToColor, ( XtPointer ) NULL);

            oofeg_add_button("TOGGLE_SCALE_GREY", "Grey scale",
                             commandWidgetClass, color_scale_palette, setIcsToGrey, ( XtPointer ) NULL);

            oofeg_add_palette("< Setup Grey Scale Colors>", color_scale_palette, & grey_scale_setup_palette);
            {
                /* CUSTOM/VARIABLE_PLOT/COLOR_SCALE_SETUP */
                Widget lx, ly;
                int ac;
                Arg al [ 6 ];

                grey_scale_setup_palette = XtCreateManagedWidget("grey_scale_setup_form", formWidgetClass,
                                                                 grey_scale_setup_palette, NULL, 0);
                ac = 0;
                XtSetArg(al [ ac ], XtNlabel, " grey_min (>=0.0)");
                ac++;
                lx = XtCreateManagedWidget("smin", labelWidgetClass, grey_scale_setup_palette, al, ac);
                ac = 0;
                XtSetArg(al [ ac ], ( String ) XtNfromHoriz, lx);
                ac++;
                XtSetArg(al [ ac ], XtNeditType, XawtextEdit);
                ac++;
                greyscale_min = oofeg_add_button("greyscale_min_val", "",
                                                 asciiTextWidgetClass, grey_scale_setup_palette,
                                                 al, ac, NULL, NULL);
                ac = 0;
                XtSetArg(al [ ac ], ( String ) XtNfromHoriz, greyscale_min);
                ac++;
                XtSetArg(al [ ac ], XtNlabel, " grey_max (<=1.0)");
                ac++;
                ly = XtCreateManagedWidget("smax", labelWidgetClass, grey_scale_setup_palette, al, ac);
                ac = 0;
                XtSetArg(al [ ac ], ( String ) XtNfromHoriz, ly);
                ac++;
                XtSetArg(al [ ac ], XtNeditType, XawtextEdit);
                ac++;
                greyscale_max = oofeg_add_button("greyscale_max_val", "",
                                                 asciiTextWidgetClass, grey_scale_setup_palette,
                                                 al, ac, NULL, NULL);

                ac = 0;
                XtSetArg(al [ ac ], ( String ) XtNfromHoriz, greyscale_max);
                ac++;
                scale_setup_ok = oofeg_add_button("precinput_xyz_ok", " OK ",
                                                  commandWidgetClass, grey_scale_setup_palette,
                                                  al, ac, pass_setgreyscale_command, ( XtPointer ) 0);
                XtAppAddActions( XtWidgetToApplicationContext( ESITopLevelWidget() ),
                                oofeg_remap_return, XtNumber(oofeg_remap_return) );
                tt1 = XtParseTranslationTable("#override <KeyPress>Return: oofegretActCmd(0)");
                XtOverrideTranslations(greyscale_min, tt1);
                XtOverrideTranslations(greyscale_max, tt1);
            }

            oofeg_add_button("TOGGLE_SMOOTH_SCALE", "Toggle smooth colors (HS)",
                             commandWidgetClass, color_scale_palette, toggleSmoothScale, ( XtPointer ) NULL);

            oofeg_add_button("TOGGLE_TRANSPARENT_CONTOURS", "Toggle transparent contours (HS)",
                             commandWidgetClass, color_scale_palette, toggleTransparentContours, ( XtPointer ) NULL);

            oofeg_add_dialog(NULL, "Set Contour Width", "Set contour width", "0", color_scale_palette,
                             apply_change, "CONTOUR_WIDTH", ESIDialogValueNumber, NULL);


            oofeg_add_dialog(NULL, "Set number of contours",
                             "Set number of contours (requires Hidden or Shaded mode)", "12",
                             color_scale_palette,
                             apply_change, "NUMBER_OF_CONTOURS", ESIDialogValueNumber, NULL);



            oofeg_add_button("TOGGLE_AUTOSCALE", "Autoscale",
                             commandWidgetClass, color_scale_palette, setAutoscale, ( XtPointer ) NULL);


            oofeg_add_palette("< Setup Color Scale >", color_scale_palette, & color_scale_setup_palette);
            {
                /* CUSTOM/VARIABLE_PLOT/COLOR_SCALE_SETUP */

                Widget lx, ly;
                int ac;
                Arg al [ 6 ];

                color_scale_setup_palette = XtCreateManagedWidget("color_scale_setup_form", formWidgetClass,
                                                                  color_scale_setup_palette, NULL, 0);
                ac = 0;
                XtSetArg(al [ ac ], XtNlabel, " scale_min");
                ac++;
                lx = XtCreateManagedWidget("smin", labelWidgetClass, color_scale_setup_palette, al, ac);
                ac = 0;
                XtSetArg(al [ ac ], ( String ) XtNfromHoriz, lx);
                ac++;
                XtSetArg(al [ ac ], XtNeditType, XawtextEdit);
                ac++;
                scale_min = oofeg_add_button("scale_min_val", "",
                                             asciiTextWidgetClass, color_scale_setup_palette,
                                             al, ac, NULL, NULL);
                ac = 0;
                XtSetArg(al [ ac ], ( String ) XtNfromHoriz, scale_min);
                ac++;
                XtSetArg(al [ ac ], XtNlabel, " scale_max");
                ac++;
                ly = XtCreateManagedWidget("smax", labelWidgetClass, color_scale_setup_palette, al, ac);
                ac = 0;
                XtSetArg(al [ ac ], ( String ) XtNfromHoriz, ly);
                ac++;
                XtSetArg(al [ ac ], XtNeditType, XawtextEdit);
                ac++;
                scale_max = oofeg_add_button("scale_max_val", "",
                                             asciiTextWidgetClass, color_scale_setup_palette,
                                             al, ac, NULL, NULL);

                ac = 0;
                XtSetArg(al [ ac ], ( String ) XtNfromHoriz, scale_max);
                ac++;
                scale_setup_ok = oofeg_add_button("precinput_xyz_ok", " OK ",
                                                  commandWidgetClass, color_scale_setup_palette,
                                                  al, ac, pass_setscale_command, ( XtPointer ) 0);
                XtAppAddActions( XtWidgetToApplicationContext( ESITopLevelWidget() ),
                                oofeg_remap_return, XtNumber(oofeg_remap_return) );
                tt1 = XtParseTranslationTable("#override <KeyPress>Return: oofegretActCmd(0)");
                XtOverrideTranslations(scale_min, tt1);
                XtOverrideTranslations(scale_max, tt1);
            }
        }

        /* CUSTOM / status */
        oofeg_add_button("STATUS_TOGGLE", "Status on/off", commandWidgetClass, view_palette, toggleStatus, NULL);



        /* CUSTOM / ACTIVE_PROBLEM */
        if ( problem->giveNumberOfSlaveProblems() ) {
            oofeg_add_dialog(NULL, "Active problem", "Problem (1 <= problem  < MAX_PROB)", "1", parent_pane,
                             apply_change, "ACTIVE_PROB", ESIDialogValueNumber, NULL);
        }


        /* CUSTOM / ACTIVE_STEP */
        oofeg_add_palette("< Active step >", parent_pane, & activeStep_palette);
        oofeg_add_palette("< Layer Update >", activeStep_palette, & layers_update_palette);
        {
            for ( i = 0; i < OOFEG_LAST_LAYER; i++ ) {
                sprintf(ltname, "%s", OOFEG_layer_names [ i ]);
                sprintf(tname, "layer_toggle_%d", i);
                oofeg_add_button(tname, ltname, toggleWidgetClass, layers_update_palette, NULL, NULL);
            }

            oofeg_add_button("layers_apply", "Apply",
                             commandWidgetClass, layers_update_palette, apply_layer_update, ( XtPointer ) layers_update_palette);
        }
        oofeg_add_dialog(NULL, "SeekStep", "Step (0 <= isptep < MAX_STEPS)", "0", activeStep_palette,
                         apply_change, "ACTIVE_STEP", ESIDialogValueNumber, NULL);
        oofeg_add_button("NEXT_STEP", "NextStep", commandWidgetClass, activeStep_palette, nextStep, NULL);
        oofeg_add_button("NEXT_STEP", "PrevStep", commandWidgetClass, activeStep_palette, previousStep, NULL);
        oofeg_add_dialog(NULL, "Seek EigMode", "Set Active Eigen value & vector (0 <= index < nroot)", "0", activeStep_palette,
                         apply_change, "ACTIVE_EIGEN_VALUE", ESIDialogValueNumber, NULL);

        oofeg_add_palette("< Animate/Loop steps >", activeStep_palette, & animate_setup_palette);
        {
            /* CUSTOM/ACTIVE_STEP/Animate */

            Widget lx, ly;
            int ac;
            Arg al [ 6 ];

            animate_setup_palette = XtCreateManagedWidget("animate_setup_form", formWidgetClass,
                                                          animate_setup_palette, NULL, 0);
            ac = 0;
            XtSetArg(al [ ac ], XtNlabel, "starting step");
            ac++;
            lx = XtCreateManagedWidget("sstep", labelWidgetClass, animate_setup_palette, al, ac);
            ac = 0;
            XtSetArg(al [ ac ], ( String ) XtNfromHoriz, lx);
            ac++;
            XtSetArg(al [ ac ], XtNeditType, XawtextEdit);
            ac++;
            start_step = oofeg_add_button("start_step", "",
                                          asciiTextWidgetClass, animate_setup_palette,
                                          al, ac, NULL, NULL);
            ac = 0;
            XtSetArg(al [ ac ], ( String ) XtNfromHoriz, start_step);
            ac++;
            XtSetArg(al [ ac ], XtNlabel, " end step");
            ac++;
            ly = XtCreateManagedWidget("estep", labelWidgetClass, animate_setup_palette, al, ac);
            ac = 0;
            XtSetArg(al [ ac ], ( String ) XtNfromHoriz, ly);
            ac++;
            XtSetArg(al [ ac ], XtNeditType, XawtextEdit);
            ac++;
            end_step = oofeg_add_button("end_step", "",
                                        asciiTextWidgetClass, animate_setup_palette,
                                        al, ac, NULL, NULL);

            ac = 0;
            XtSetArg(al [ ac ], ( String ) XtNfromHoriz, end_step);
            ac++;
            animate_scale_setup_ok = oofeg_add_button("precinput_xyz_ok", " OK ",
                                                      commandWidgetClass, animate_setup_palette,
                                                      al, ac, pass_setanimate_command, ( XtPointer ) 0);
            XtAppAddActions( XtWidgetToApplicationContext( ESITopLevelWidget() ),
                            oofeg_remap_return, XtNumber(oofeg_remap_return) );
            tt1 = XtParseTranslationTable("#override <KeyPress>Return: oofegretActCmd(0)");
            XtOverrideTranslations(scale_min, tt1);
            XtOverrideTranslations(scale_max, tt1);
        }

        /* CUSTOM / PLOT */
        oofeg_add_palette("< Mesh Plot >", parent_pane, & plot_palette);
        oofeg_add_button("GEO_PLOT", "Geom Plot", commandWidgetClass, plot_palette, GeoPlot, NULL);

        oofeg_add_palette("< DofMan Plot >", plot_palette, & dofman_menu);
        oofeg_add_button("DOFMAN_PLOT", "Node Geometry", commandWidgetClass, dofman_menu, nodePlot, NULL);
        oofeg_add_button("DOFMANNUM_PLOT", "DofMan Numbers", commandWidgetClass, dofman_menu, nodeAnnPlot, NULL);
        oofeg_add_button("ELEMNUM_PLOT", "Element Numbers", commandWidgetClass, dofman_menu, elementAnnPlot, NULL);

        oofeg_add_button( "BCE_PLOT", "essential BC", commandWidgetClass, dofman_menu, bcPlot, ( XtPointer ) ( vectorAddr + 0 ) );
        oofeg_add_button( "BCN_PLOT", "natural   BC", commandWidgetClass, dofman_menu, bcPlot, ( XtPointer ) ( vectorAddr + 1 ) );

        oofeg_add_palette("< Deformations >", plot_palette, & deformations_menu);
        oofeg_add_button("DEF_PLOT", "DefGeom Plot", commandWidgetClass, deformations_menu, defPlot, NULL);
        oofeg_add_button("EIGVEC_PLOT", "EigVec Plot", commandWidgetClass, deformations_menu, eigVecPlot, NULL);
        oofeg_add_dialog(NULL, "Set DefScale", "Input Deformation Scale (0 <= scale < MAX_Scale)", "0.0", deformations_menu,
                         apply_change, "SET_DEF_SCALE", ESIDialogValueNumber, NULL);
        oofeg_add_button("AutoSetDefScale", "Auto DefScale", commandWidgetClass, deformations_menu, defAutoScale, NULL);
        oofeg_add_button("SPARSE_PLOT", "StiffSparse Plot", commandWidgetClass, plot_palette,
                         showSparseMtrxStructure, NULL);

        /* CUSTOM/ internal vars plot */
        oofeg_add_palette("< Variable plot >", parent_pane, & varplot_palette);
        {
            oofeg_add_button("DEFORMED_PLOT_FLAG", "Plot on Deformed shape", toggleWidgetClass, varplot_palette, NULL, NULL);

            oofeg_add_popdown_menu("< Data Mode >", varplot_palette, & valmode_palette);
            {
                oofeg_add_menu_item( "VALUE_MODE_RECOVERED", "Recovered Values", valmode_palette, valmodesel, ( XtPointer ) ( vectorAddr + 0 ) );
                oofeg_add_menu_item( "VALUE_MODE_LOCAL", "Local Values", valmode_palette, valmodesel, ( XtPointer ) ( vectorAddr + 1 ) );
            }

            oofeg_add_popdown_menu("< Plot Algorithm >", varplot_palette, & plotalgo_palette);
            {
                oofeg_add_menu_item( "ISO_SURF_PLOT", "IsoSurf Plot", plotalgo_palette, plotAlgosel, ( XtPointer ) ( vectorAddr + 0 ) );
                oofeg_add_menu_item( "ISO_LINE_PLOT", "IsoLine Plot", plotalgo_palette, plotAlgosel, ( XtPointer ) ( vectorAddr + 1 ) );
                oofeg_add_menu_item( "ISO_LANDPROFILE", "Z-Profile Plot", plotalgo_palette, plotAlgosel, ( XtPointer ) ( vectorAddr + 2 ) );
                oofeg_add_menu_item( "ISO_LANDCOLORPROFILE", "Z-ColorProfile Plot", plotalgo_palette,
                                    plotAlgosel, ( XtPointer ) ( vectorAddr + 3 ) );
            }

            oofeg_add_palette("< Scalar plot >", varplot_palette, & scalarplot_palette);
            {
                oofeg_add_button("CYLINDRICAL_CS_FLAG", "Cylindrical cs", toggleWidgetClass, scalarplot_palette, NULL, NULL);
                oofeg_add_popdown_menu("< Stress/Force plot >", scalarplot_palette, & stressForce_palette);
                {
                    oofeg_add_menu_item( "SX_PLOT", "Sxx Stress/Force Plot", stressForce_palette, stresscompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "SY_PLOT", "Syy Stress/Force Plot", stressForce_palette, stresscompPlot, ( XtPointer ) ( vectorAddr + 1 ) );
                    oofeg_add_menu_item( "SZ_PLOT", "Szz Stress/Force Plot", stressForce_palette, stresscompPlot, ( XtPointer ) ( vectorAddr + 2 ) );
                    oofeg_add_menu_item( "SYZ_PLOT", "Syz Stress/Force Plot", stressForce_palette, stresscompPlot, ( XtPointer ) ( vectorAddr + 3 ) );
                    oofeg_add_menu_item( "SZX_PLOT", "Szx Stress/Force Plot", stressForce_palette, stresscompPlot, ( XtPointer ) ( vectorAddr + 4 ) );
                    oofeg_add_menu_item( "SXY_PLOT", "Sxy Stress/Force Plot", stressForce_palette, stresscompPlot, ( XtPointer ) ( vectorAddr + 5 ) );
                    oofeg_add_menu_item("SEPARATOR", "--------------------", stressForce_palette, emptyaction, NULL);
                    oofeg_add_menu_item( "S11_PLOT", "S11 Stress/Force Plot", stressForce_palette,
                                        princstresscompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "S22_PLOT", "S22 Stress/Force Plot", stressForce_palette,
                                        princstresscompPlot, ( XtPointer ) ( vectorAddr + 1 ) );
                    oofeg_add_menu_item( "S33_PLOT", "S33 Stress/Force Plot", stressForce_palette,
                                        princstresscompPlot, ( XtPointer ) ( vectorAddr + 2 ) );
                }
                oofeg_add_popdown_menu("< Strain plot >", scalarplot_palette, & scalarstrain_palette);
                {
                    oofeg_add_menu_item( "SX_PLOT", "Eps_x ", scalarstrain_palette, straincompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "SY_PLOT", "Eps_y", scalarstrain_palette, straincompPlot, ( XtPointer ) ( vectorAddr + 1 ) );
                    oofeg_add_menu_item( "SZ_PLOT", "Eps_z", scalarstrain_palette, straincompPlot, ( XtPointer ) ( vectorAddr + 2 ) );
                    oofeg_add_menu_item( "SYZ_PLOT", "Gam_yz", scalarstrain_palette, straincompPlot, ( XtPointer ) ( vectorAddr + 3 ) );
                    oofeg_add_menu_item( "SZX_PLOT", "Gam_xz", scalarstrain_palette, straincompPlot, ( XtPointer ) ( vectorAddr + 4 ) );
                    oofeg_add_menu_item( "SXY_PLOT", "Gam_xy", scalarstrain_palette, straincompPlot, ( XtPointer ) ( vectorAddr + 5 ) );
                    oofeg_add_menu_item("SEPARATOR", "--------------------", scalarstrain_palette, emptyaction, NULL);
                    oofeg_add_menu_item( "E11_PLOT", "Eps_11 ", scalarstrain_palette, princstraincompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "E22_PLOT", "Eps_22", scalarstrain_palette, princstraincompPlot, ( XtPointer ) ( vectorAddr + 1 ) );
                    oofeg_add_menu_item( "E33_PLOT", "Eps_33", scalarstrain_palette, princstraincompPlot, ( XtPointer ) ( vectorAddr + 2 ) );
                    oofeg_add_menu_item("SEPARATOR", "--------------------", scalarstrain_palette, emptyaction, NULL);
                    oofeg_add_menu_item( "EES_PLOT", "Equiv_Strain", scalarstrain_palette, equivstraincompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                }

                oofeg_add_popdown_menu("< Plastic strain plot >", scalarplot_palette, & scalarplasticstrain_palette);
                {
                    oofeg_add_menu_item( "SX_PLOT", "Eps_x ", scalarplasticstrain_palette, plaststraincompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "SY_PLOT", "Eps_y", scalarplasticstrain_palette, plaststraincompPlot, ( XtPointer ) ( vectorAddr + 1 ) );
                    oofeg_add_menu_item( "SZ_PLOT", "Eps_z", scalarplasticstrain_palette, plaststraincompPlot, ( XtPointer ) ( vectorAddr + 2 ) );
                    oofeg_add_menu_item( "SYZ_PLOT", "Gam_yz", scalarplasticstrain_palette, plaststraincompPlot, ( XtPointer ) ( vectorAddr + 3 ) );
                    oofeg_add_menu_item( "SZX_PLOT", "Gam_xz", scalarplasticstrain_palette, plaststraincompPlot, ( XtPointer ) ( vectorAddr + 4 ) );
                    oofeg_add_menu_item( "SXY_PLOT", "Gam_xy", scalarplasticstrain_palette, plaststraincompPlot, ( XtPointer ) ( vectorAddr + 5 ) );
                    oofeg_add_menu_item("SEPARATOR", "--------------------", scalarplasticstrain_palette, emptyaction, NULL);
                    oofeg_add_menu_item( "E11_PLOT", "Eps_11 ", scalarplasticstrain_palette,
                                        princplaststraincompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "E22_PLOT", "Eps_22", scalarplasticstrain_palette,
                                        princplaststraincompPlot, ( XtPointer ) ( vectorAddr + 1 ) );
                    oofeg_add_menu_item( "E33_PLOT", "Eps_33", scalarplasticstrain_palette,
                                        princplaststraincompPlot, ( XtPointer ) ( vectorAddr + 2 ) );
                }

                oofeg_add_popdown_menu("< Beam Forces/Moments>", scalarplot_palette, & bendingMoment_palette);
                {
                    oofeg_add_menu_item( "NX_PLOT", "Nx Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "NY_PLOT", "Ny Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 1 ) );
                    oofeg_add_menu_item( "NZ_PLOT", "Nz Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 2 ) );
                    oofeg_add_menu_item( "VYZ_PLOT", "Vyz Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 3 ) );
                    oofeg_add_menu_item( "VXZ_PLOT", "Vxz Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 4 ) );
                    oofeg_add_menu_item( "VXY_PLOT", "Vxy Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 5 ) );


                    oofeg_add_menu_item( "MX_PLOT", "Mx Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 6 ) );
                    oofeg_add_menu_item( "MY_PLOT", "My Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 7 ) );
                    oofeg_add_menu_item( "MZ_PLOT", "Mz Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 8 ) );
                    oofeg_add_menu_item( "MYZ_PLOT", "Myz Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 9 ) );
                    oofeg_add_menu_item( "MXZ_PLOT", "Mxz Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 10 ) );
                    oofeg_add_menu_item( "MXY_PLOT", "Mxy Plot", bendingMoment_palette, beamForcePlot, ( XtPointer ) ( vectorAddr + 11 ) );
                }

                oofeg_add_popdown_menu("< Plot Damage >", scalarplot_palette, & scalardamage_palette);
                {
                    oofeg_add_menu_item( "DX_PLOT", "D_x ", scalardamage_palette, damagecompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "DY_PLOT", "D_y", scalardamage_palette, damagecompPlot, ( XtPointer ) ( vectorAddr + 1 ) );
                    oofeg_add_menu_item( "DZ_PLOT", "D_z", scalardamage_palette, damagecompPlot, ( XtPointer ) ( vectorAddr + 2 ) );
                    oofeg_add_menu_item( "DYZ_PLOT", "D_yz", scalardamage_palette, damagecompPlot, ( XtPointer ) ( vectorAddr + 3 ) );
                    oofeg_add_menu_item( "DZX_PLOT", "D_zz", scalardamage_palette, damagecompPlot, ( XtPointer ) ( vectorAddr + 4 ) );
                    oofeg_add_menu_item( "DXY_PLOT", "D_xy", scalardamage_palette, damagecompPlot, ( XtPointer ) ( vectorAddr + 5 ) );
                    oofeg_add_menu_item("SEPARATOR", "--------------------", scalardamage_palette, emptyaction, NULL);
                    oofeg_add_menu_item( "D11_PLOT", "D_11 ", scalardamage_palette, princdamagecompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "D22_PLOT", "D_22", scalardamage_palette, princdamagecompPlot, ( XtPointer ) ( vectorAddr + 1 ) );
                    oofeg_add_menu_item( "D33_PLOT", "D_33", scalardamage_palette, princdamagecompPlot, ( XtPointer ) ( vectorAddr + 2 ) );
                }

                //mj
                oofeg_add_popdown_menu("< Internal Vars >", scalarplot_palette, & scalardamage_palette);
                {
                    oofeg_add_menu_item( "EQUIV_STRAIN", "eps_eq ", scalardamage_palette, epseqcompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "CUM_PLAST_STRAIN", "kappa ", scalardamage_palette, kappacompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "CUM_PLAST_STRAIN_2", "kappa_nl ", scalardamage_palette, kappa2compPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "DISS_WORK", "diss_work ", scalardamage_palette, dissworkcompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "STRESS_WORK", "stress_work ", scalardamage_palette, stressworkcompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "FREE_ENERGY", "free_energy ", scalardamage_palette, freeenergycompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                }
                //emj

                oofeg_add_button( "TEMPERATURE_PLOT", "Temperature field", commandWidgetClass, scalarplot_palette,
                                 varPlot, ( XtPointer ) ( oofeg_draw_modes + 17 ) );

                oofeg_add_button( "CONCENTRATION1_PLOT", "MassConcentration(1)", commandWidgetClass, scalarplot_palette,
                                 varPlot, ( XtPointer ) ( oofeg_draw_modes + 18 ) );

                oofeg_add_button( "PRESSURE_PLOT", "Pressure", commandWidgetClass, scalarplot_palette,
                                 varPlot, ( XtPointer ) ( oofeg_draw_modes + 20 ) );

                oofeg_add_button( "VOF_PLOT", "VOF Fraction", commandWidgetClass, scalarplot_palette,
                                 varPlot, ( XtPointer ) ( oofeg_draw_modes + 21 ) );

                oofeg_add_button( "DENSITY_PLOT", "Density", commandWidgetClass, scalarplot_palette,
                                 varPlot, ( XtPointer ) ( oofeg_draw_modes + 22 ) );


                oofeg_add_popdown_menu("< Error plot >", scalarplot_palette, & scalarerror_palette);
                {
                    oofeg_add_menu_item( "Error_Indicator", "Error Indicator ", scalarerror_palette,
                                        errorcompPlot, ( XtPointer ) ( vectorAddr + 0 ) );
                    oofeg_add_menu_item( "ERROR_STRESSINTERRIOR", "Stress Error", scalarerror_palette,
                                        errorcompPlot, ( XtPointer ) ( vectorAddr + 1 ) );
                    oofeg_add_menu_item( "ERROR_PRIMARYUNKNOWN", "Unknown Error", scalarerror_palette,
                                        errorcompPlot, ( XtPointer ) ( vectorAddr + 2 ) );
                }

                oofeg_add_button( "MESHDENS_PLOT", "RelMesh Density", commandWidgetClass, scalarplot_palette,
                                 varPlot, ( XtPointer ) ( oofeg_draw_modes + 16 ) );

                oofeg_add_dialog(NULL, "Set ZProfScale", "Input z-Profile Scale (0 <= scale < MAX_Scale)", "0.0", scalarplot_palette,
                                 apply_change, "SET_ZPROF_SCALE", ESIDialogValueNumber, NULL);
            }

            oofeg_add_palette("< Vector plot >", varplot_palette, & vectorplot_palette);
            oofeg_add_button( "VELOCITY_PLOT", "Velocity", commandWidgetClass, vectorplot_palette,
                             varPlot, ( XtPointer ) ( oofeg_draw_modes + 19 ) );

            oofeg_add_palette("< Tensor plot >", varplot_palette, & tensorplot_palette);


            oofeg_add_button( "YIELD_PLOT", "Plot YieldedElems",
                             commandWidgetClass, varplot_palette, varPlot, ( XtPointer ) ( oofeg_draw_modes + 12 ) );

            oofeg_add_button( "CRACK_PLOT", "Plot CrackPattern",
                             commandWidgetClass, varplot_palette, varPlot, ( XtPointer ) ( oofeg_draw_modes + 13 ) );
        }

        oofeg_add_palette("< Set smoother >", varplot_palette, & smoother_palette);
        {
            oofeg_add_button( "NodalAveraging", "NodalAveraging Smoother",
                             commandWidgetClass, smoother_palette, setSmoother, ( XtPointer ) ( oofeg_smoother_modes + 0 ) );

            oofeg_add_button( "ZieZhu", "ZZ Smoother",
                             commandWidgetClass, smoother_palette, setSmoother, ( XtPointer ) ( oofeg_smoother_modes + 1 ) );

            oofeg_add_button( "SPR", "SPR Smoother",
                             commandWidgetClass, smoother_palette, setSmoother, ( XtPointer ) ( oofeg_smoother_modes + 2 ) );
        }
    }

    /* CUSTOM/FILTERS */
    oofeg_add_palette("< Filters >", parent_pane, & filters_palette);
    oofeg_add_palette("< Material Region Filter >", filters_palette, & matregfilter_palette);
    for ( int id = 1; id <= problem->giveNumberOfDomains(); id++ ) {
        nmat = problem->giveDomain(id)->giveNumberOfMaterialModels();

        for ( i = 1; i <= nmat; i++ ) {
            sprintf(ltname, "region %2d.%2d", id, i);
            sprintf(tname, "mat_reg_toggle_%d-%d", id, i);
            oofeg_add_button(tname, ltname, toggleWidgetClass, matregfilter_palette, NULL, NULL);
        }
    }

    oofeg_add_button("mat_reg_filter_apply", "Apply",
                     commandWidgetClass, matregfilter_palette, apply_mat_reg_filter, ( XtPointer ) matregfilter_palette);

    oofeg_add_dialog(NULL, "Element Filter", "Enter active element range {(s1 e1) (s2 e2) i1 i2}", "", filters_palette,
                     apply_change, "element_filter", ESIDialogValueString, NULL);


    /* CUSTOM/LAYERS */
    oofeg_add_palette("< Layers >", parent_pane, & layers_palette);
    for ( i = 0; i < OOFEG_LAST_LAYER; i++ ) {
        sprintf(ltname, "%s", OOFEG_layer_names [ i ]);
        sprintf(tname, "layer_toggle_%d", i);
        oofeg_add_button(tname, ltname, toggleWidgetClass, layers_palette, NULL, NULL);
    }

    oofeg_add_button("layers_apply", "Apply",
                     commandWidgetClass, layers_palette, apply_layer_settings, ( XtPointer ) layers_palette);


    /* CUSTOM / info */
    oofeg_add_button("INFO", "Info", commandWidgetClass, parent_pane, print_state, NULL);

#ifdef OOFEG_DEVEL
    /* CUSTOM / DEBUG */
    oofeg_add_palette("< DEBUG >", parent_pane, & debug_palette);
    {
        oofeg_add_button("Run job", "Run job",
                         commandWidgetClass, debug_palette, debug_run, NULL);
    }
#endif
}

static void apply_change(Widget wid, XtPointer cl, XtPointer cd)
{
    char buf [ 512 ];

    if ( cd != NULL ) {
        sprintf(buf, "%s %s", ( char * ) cl, ( char * ) cd);
    } else {
        sprintf(buf, "%s", ( char * ) cl);
    }

    OOFEGSimpleCmd(buf);
}

void OOFEGSimpleCmd(char *buf)
{
    char cmd [ 20 ];
    char *remain;
    int stepinfo [ 2 ];
    int istep, pstep, iversion = 0;
    contextIOResultType result;

    readSimpleString(buf, cmd, & remain); // read comand
    if ( !strncasecmp(cmd, "active_step", 11) ) {
        pstep = gc [ 0 ].getActiveStep();
        istep = atoi(remain);
        stepinfo [ 0 ] = istep;
        stepinfo [ 1 ] = iversion;
        try {
            result = problem->restoreContext(NULL, CM_State, ( void * ) stepinfo);
        } catch(ContextIOERR & m) {
            m.print();
            stepinfo [ 0 ] = pstep;
            stepinfo [ 1 ] = iversion;
            try {
                result = problem->restoreContext(NULL, CM_State, ( void * ) stepinfo);
            } catch(ContextIOERR & m2) {
                m2.print();
                exit(1);
            }
            istep = pstep;
        }
        gc [ 0 ].setActiveStep(istep);

        updateGraphics();

        //problem ->forceEquationNumbering();
    } else if ( !strncasecmp(cmd, "active_eigen_value", 18) ) {
        istep = atoi(remain);
        gc [ 0 ].setActiveStep(istep);
        stepinfo [ 0 ] = istep;
        stepinfo [ 1 ] = iversion;
        try {
            problem->restoreContext(NULL, CM_State, ( void * ) stepinfo);
        } catch(ContextIOERR & m) {
            m.print();
            exit(1);
        }
        gc [ 0 ].setActiveEigVal(istep);
    } else if ( !strncasecmp(cmd, "set_def_scale", 13) ) {
        gc [ 0 ].setDefScale( strtod(remain, NULL) );
        if ( strtod(remain, NULL) < 0 ) {
            gc [ 0 ].setDefScale(1.0);
        }
    } else if ( !strncasecmp(cmd, "set_zprof_scale", 15) ) {
        gc [ 0 ].setLandScale( strtod(remain, NULL) );
        if ( strtod(remain, NULL) < 0 ) {
            gc [ 0 ].setDefScale(1.0);
        }
    } else if ( !strncasecmp(cmd, "element_filter", 14) ) {
        char buff [ 80 ];
        sprintf(buff, "element_filter %s", remain);
        gc [ 0 ].setElementFilterState(buff);
    } else if ( !strncasecmp(cmd, "active_prob", 11) ) {
        int iprob = atoi(remain);
        if ( ( iprob > 0 ) && ( iprob <= problem->giveNumberOfSlaveProblems() ) ) {
            gc [ 0 ].setActiveProblem(iprob);
            printf("Setting problem %d as active\n", iprob);
            updateGraphics();
        }
    } else if ( !strncasecmp(cmd, "contour_width", 13) ) {
        int width = atoi(remain);
        SetContourWidth(width);

        EView *v = ( EView * ) get_list_next(age_model->dependent_views, NULL);
        while ( v != NULL ) {
            EVFastRedraw(v);
            v = ( EView * ) get_list_next(age_model->dependent_views, v);
        }
    } else if ( !strncasecmp(cmd, "number_of_contours", 18) ) {
        SetContourCount( atoi(remain) );
    }
}


void GeoPlot(Widget wid, XtPointer cl, XtPointer cd)
{
    gc [ OOFEG_RAW_GEOMETRY_LAYER ].setPlotMode(OGC_rawGeometry);
    drawData(gc [ OOFEG_RAW_GEOMETRY_LAYER ]);
}


void defPlot(Widget wid, XtPointer cl, XtPointer cd)
{
    if ( gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].getActiveStep() != -1 ) {
        gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].setPlotMode(OGC_deformedGeometry);
        drawData(gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ]);
    }
}


void  defAutoScale(Widget wid, XtPointer cl, XtPointer cd)
{
    int init = 1, id, i, j, nnodes;
    double mincoords [ 3 ], maxcoords [ 3 ];
    double coord, maxdef = 0.0;
    Domain *domain;
    DofManager *dman;

    TimeStep *tStep = gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].getActiveProblem()->giveCurrentStep();
    if ( tStep == NULL ) {
        gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].setDefScale(0.0);
        return;
    }

    for ( id = 1; id <= gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].getActiveProblem()->giveNumberOfDomains(); id++ ) {
        domain = gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].getActiveProblem()->giveDomain(id);
        nnodes = domain->giveNumberOfDofManagers();
        if ( id == 1 ) {
            for ( j = 1; j <= 3; j++ ) {
                maxcoords [ j - 1 ] = mincoords [ j - 1 ] = domain->giveNode(1)->giveCoordinate(j);
            }
        }

        for ( i = 1; i <= nnodes; i++ ) {
            dman = domain->giveDofManager(i);
            if ( ( dman->giveClassID() == NodeClass ) || ( dman->giveClassID() == RigidArmNodeClass ) ) {
                if ( init ) {
                    for ( j = 1; j <= 3; j++ ) {
                        maxcoords [ j - 1 ] = mincoords [ j - 1 ] = domain->giveNode(i)->giveCoordinate(j);
                    }

                    init = 0;
                }

                for ( j = 1; j <= 3; j++ ) {
                    coord = domain->giveNode(i)->giveCoordinate(j);
                    maxdef = max( maxdef, fabs(domain->giveNode(i)->giveUpdatedCoordinate(j, tStep, EID_MomentumBalance, 1.0) - coord) );
                    maxcoords [ j - 1 ] = max(maxcoords [ j - 1 ], coord);
                    mincoords [ j - 1 ] = min(mincoords [ j - 1 ], coord);
                }
            }
        }
    }

    double maxCoordDiff = maxcoords [ 0 ] - mincoords [ 0 ];
    maxCoordDiff = max(maxCoordDiff, maxcoords [ 1 ] - mincoords [ 1 ]);
    maxCoordDiff = max(maxCoordDiff, maxcoords [ 2 ] - mincoords [ 2 ]);

    gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].setDefScale(0.1 * maxCoordDiff / maxdef);
    printf("Def scale %e\n", 0.1 * maxCoordDiff / maxdef);
}




void eigVecPlot(Widget wid, XtPointer cl, XtPointer cd)
{
    if ( gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].getActiveStep() != -1 ) {
        gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].setPlotMode(OGC_eigenVectorGeometry);
        drawData(gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ]);
    }
}


int  updateDefPlotFlag()
{
    char tname [ 32 ];
    Widget t;
    Arg al [ 1 ];
    int ac;
    Boolean s;

    sprintf(tname, "DEFORMED_PLOT_FLAG");
    t = XtNameToWidget(varplot_palette, tname);
    if ( t != NULL ) {
        ac = 0;
        XtSetArg(al [ ac ], ( String ) XtNstate, & s);
        ac++;
        XtGetValues(t, al, ac);
        gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].setInternalVarsDefGeoFlag( ( BOOLEAN ) s );
    }

    return 1;
}





void nextStep(Widget wid, XtPointer cl, XtPointer cd)
{
    contextIOResultType result;
    int istep, prevStep, stepStep, prevStepVersion, istepVersion;
    int stepInfo [ 2 ];

    stepStep = problem->giveContextOutputStep();
    if ( stepStep == 0 ) {
        stepStep = 1;
    }

    if ( gc [ 0 ].getActiveStep() != -1 ) {
        prevStep = gc [ 0 ].getActiveStep();
        prevStepVersion = gc [ 0 ].getActiveStepVersion();

        // first try next version for the same step
        if ( problem->testContextFile(prevStep, prevStepVersion + 1) ) {
            istepVersion = prevStepVersion + 1;
            stepInfo [ 0 ] = prevStep;
            stepInfo [ 1 ] = istepVersion;
            printf("OOFEG: restoring context file %d.%d\n", stepInfo [ 0 ], stepInfo [ 1 ]);
            try {
                result = problem->restoreContext(NULL, CM_State, ( void * ) stepInfo);
            } catch(ContextIOERR & m) {
                m.print();
                istepVersion = 0;
                stepInfo [ 0 ] = prevStep;
                stepInfo [ 1 ] = 0;
                try {
                    problem->restoreContext(NULL, CM_State, ( void * ) stepInfo);
                } catch(ContextIOERR & m2) {
                    m2.print();
                    exit(1);
                }
            }
            gc [ 0 ].setActiveStep(prevStep);
            gc [ 0 ].setActiveStepVersion(istepVersion);
        } else {
            // no next version exist => next step with version 0
            istep = prevStep + stepStep;
            stepInfo [ 0 ] = istep;
            stepInfo [ 1 ] = 0;

            //printf ("NextStep: prevStep %d, nstep %d, stepStep %d\n", prevStep, istep, stepStep);
            try {
                result = problem->restoreContext(NULL, CM_State, ( void * ) stepInfo);
            } catch(ContextIOERR & m) {
                m.print();
                stepInfo [ 0 ] = prevStep;
                stepInfo [ 1 ] = 0;
                try {
                    problem->restoreContext(NULL, CM_State, ( void * ) stepInfo);
                } catch(ContextIOERR & m2) {
                    m2.print();
                    exit(1);
                }
                istep = prevStep;
            }
            gc [ 0 ].setActiveStep(istep);
            gc [ 0 ].setActiveStepVersion(0);
        }
    } else {
        istep = problem->giveNumberOfFirstStep() + stepStep - 1;
        gc [ 0 ].setActiveStep(istep);
        stepInfo [ 0 ] = istep;
        stepInfo [ 1 ] = 0;
        try {
            result = problem->restoreContext(NULL, CM_State, ( void * ) stepInfo);
        } catch(ContextIOERR & m) {
            m.print();
            exit(1);
        }
    }

    updateGraphics();
}

void previousStep(Widget wid, XtPointer cl, XtPointer cd)
{
    contextIOResultType result;
    int istep, prevStep, stepStep = problem->giveContextOutputStep();
    int stepInfo [ 2 ];
    if ( stepStep == 0 ) {
        stepStep = 1;
    }

    if ( gc [ 0 ].getActiveStep() != -1 ) {
        prevStep = gc [ 0 ].getActiveStep();
        istep = prevStep - stepStep;
        if ( istep >= 0 ) {
            stepInfo [ 0 ] = istep;
            stepInfo [ 1 ] = 0;
            try {
                result = problem->restoreContext(NULL, CM_State, ( void * ) stepInfo);
            } catch(ContextIOERR & m) {
                m.print();
                stepInfo [ 0 ] = prevStep;
                stepInfo [ 1 ] = 0;
                try {
                    problem->restoreContext(NULL, CM_State, ( void * ) stepInfo);
                } catch(ContextIOERR & m2) {
                    m2.print();
                    exit(1);
                }
                istep = prevStep;
            }
            gc [ 0 ].setActiveStep(istep);
            gc [ 0 ].setActiveStepVersion(0);
        }
    } else {
        istep = 1;
        gc [ 0 ].setActiveStep(istep);
        gc [ 0 ].setActiveStepVersion(0);
        stepInfo [ 0 ] = istep;
        stepInfo [ 1 ] = 0;
        try {
            result = problem->restoreContext(NULL, CM_State, ( void * ) stepInfo);
        } catch(ContextIOERR & m) {
            m.print();
            exit(1);
        }
    }

    updateGraphics();
}

void nodePlot(Widget wid, XtPointer cl, XtPointer cd)
{
    gc [ OOFEG_NODE_ANNOTATION_LAYER ].setPlotMode(OGC_nodeGeometry);
    drawData(gc [ OOFEG_NODE_ANNOTATION_LAYER ]);
}

void nodeAnnPlot(Widget wid, XtPointer cl, XtPointer cd)
{
    gc [ OOFEG_NODE_ANNOTATION_LAYER ].setPlotMode(OGC_nodeAnnotation);
    drawData(gc [ OOFEG_NODE_ANNOTATION_LAYER ]);
}

void elementAnnPlot(Widget wid, XtPointer cl, XtPointer cd)
{
    gc [ OOFEG_ELEMENT_ANNOTATION_LAYER ].setPlotMode(OGC_elementAnnotation);
    drawData(gc [ OOFEG_ELEMENT_ANNOTATION_LAYER ]);
}

void bcPlot(Widget wid, XtPointer cl, XtPointer cd)
{
    int mode = * ( ( int * ) cl );
    if ( mode == 1 ) { //OGC_essentialBC
        gc [ OOFEG_BCIC_ANNOTATION_LAYER ].setPlotMode(OGC_essentialBC);
        drawData(gc [ OOFEG_BCIC_ANNOTATION_LAYER ]);
    } else if ( mode == 2 ) { // OGC_naturalBC
        gc [ OOFEG_NATURALBC_LAYER ].setPlotMode(OGC_naturalBC);
        drawData(gc [ OOFEG_NATURALBC_LAYER ]);
    }
}

void emptyaction(Widget w, XtPointer ptr, XtPointer call_data)
{ }

void stresscompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;
    Boolean s;
    char tname [ 32 ];
    Widget t;
    Arg al [ 1 ];

    updateDefPlotFlag();

    sprintf(tname, "CYLINDRICAL_CS_FLAG");
    t = XtNameToWidget(scalarplot_palette, tname);
    if ( t != NULL ) {
        ac = 0;
        XtSetArg(al [ ac ], ( String ) XtNstate, & s);
        ac++;
        XtGetValues(t, al, ac);
    }

    ac = * ( ( int * ) ptr );
    if ( s ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_CylindricalStressTensor);
    } else {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_StressTensor);
    }

    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(ac);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void princstresscompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    updateDefPlotFlag();

    ac = * ( ( int * ) ptr );
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_PrincipalStressTensor);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(ac);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void equivstraincompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    updateDefPlotFlag();

    ac = * ( ( int * ) ptr );
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_MaxEquivalentStrainLevel);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(ac);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void straincompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;
    Boolean s;
    char tname [ 32 ];
    Widget t;
    Arg al [ 1 ];

    updateDefPlotFlag();

    sprintf(tname, "CYLINDRICAL_CS_FLAG");
    t = XtNameToWidget(scalarplot_palette, tname);
    if ( t != NULL ) {
        ac = 0;
        XtSetArg(al [ ac ], ( String ) XtNstate, & s);
        ac++;
        XtGetValues(t, al, ac);
    }

    ac = * ( ( int * ) ptr );
    if ( s ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_CylindricalStrainTensor);
    } else {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_StrainTensor);
    }

    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(ac);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}


void princstraincompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    updateDefPlotFlag();

    ac = * ( ( int * ) ptr );
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_PrincipalStrainTensor);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(ac);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void plaststraincompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    updateDefPlotFlag();

    //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
    ac = * ( ( int * ) ptr );
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_PlasticStrainTensor);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(ac);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}


void princplaststraincompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    updateDefPlotFlag();

    //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
    ac = * ( ( int * ) ptr );
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_PrincipalPlasticStrainTensor);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(ac);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

//mj
void epseqcompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    updateDefPlotFlag();
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_MaxEquivalentStrainLevel);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void kappacompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    updateDefPlotFlag();
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_CumPlasticStrain);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void kappa2compPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    updateDefPlotFlag();
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_CumPlasticStrain_2);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void dissworkcompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    updateDefPlotFlag();
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_DissWorkDensity);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void stressworkcompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    updateDefPlotFlag();
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_StressWorkDensity);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void freeenergycompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    updateDefPlotFlag();
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_FreeEnergyDensity);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}
//emj

void damagecompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    updateDefPlotFlag();

    //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
    ac = * ( ( int * ) ptr );
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_DamageTensor);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(ac);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}


void princdamagecompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    updateDefPlotFlag();

    //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
    ac = * ( ( int * ) ptr );
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_PrincipalDamageTensor);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(ac);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void errorcompPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    updateDefPlotFlag();

    //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
    ac = * ( ( int * ) ptr );

    if ( ac == 1 ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_ErrorIndicatorLevel);
    } else if ( ac == 2 ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_InternalStressError);
    } else if ( ac == 3 ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_PrimaryUnknownError);
    } else {
        deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
        return;
    }

    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);
    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void beamForcePlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    updateDefPlotFlag();


    //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
    ac = * ( ( int * ) ptr );
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_ShellForceMomentumTensor);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(ac);
    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);

    setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void plotAlgosel(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    ac = * ( ( int * ) ptr );
    if ( ac == 1 ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setScalarAlgo(SA_ISO_SURF);
        SetFringeTableContourScaleOff();
    } else if ( ac == 2 ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setScalarAlgo(SA_ISO_SURF);
        SetFringeTableContourScaleOn();
    } else if ( ac == 3 ) {
        SetFringeTableContourScaleOff();
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setScalarAlgo(SA_ZPROFILE);
    } else if ( ac == 4 ) {
        SetFringeTableContourScaleOff();
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setScalarAlgo(SA_COLORZPROFILE);
    }



    gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);
    updateDefPlotFlag();
    //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void valmodesel(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;

    ac = * ( ( int * ) ptr );
    if ( ac == 1 ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateMode(ISM_recovered);
    } else if ( ac == 2 ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateMode(ISM_local);
    }

    updateDefPlotFlag();
    //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
    setupSmoother(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
        drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
    }
}

void varPlot(Widget w, XtPointer ptr, XtPointer call_data)
{
    DrawMode mode = * ( ( DrawMode * ) ptr );

    updateDefPlotFlag();

    if ( mode == damageLevel ) {
        //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);

        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_DamageTensor);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);
        setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
        if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
            drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
        }
    } else if ( mode == crackedState ) {
        //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);

        gc [ OOFEG_CRACK_PATTERN_LAYER ].setInternalStateType(IST_CrackState);
        gc [ OOFEG_CRACK_PATTERN_LAYER ].setIntVarIndx(1);
        gc [ OOFEG_CRACK_PATTERN_LAYER ].setPlotMode(OGC_elemSpecial);
        setupData(gc [ OOFEG_CRACK_PATTERN_LAYER ]);
        if ( gc [ OOFEG_CRACK_PATTERN_LAYER ].getActiveStep() != -1 ) {
            drawData(gc [ OOFEG_CRACK_PATTERN_LAYER ]);
        }

        /*
         * } else if (mode == errorIndicatorLevel) {
         *
         * deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
         * gc.setInternalStateType (IST_ErrorIndicatorLevel);
         * gc.setIntVarIndx (1);
         * gc.setPlotMode (OGC_scalarPlot);
         * setupData (gc);
         * if (gc.getActiveStep() != -1) drawData(gc);
         */
    } else if ( mode == temperatureField ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_Temperature);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);
        if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
            setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
            drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
        }
    } else if ( mode == massConcentration1Field ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_MassConcentration_1);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);
        if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
            setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
            drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
        }
    } else if ( mode == velocityField ) {
        gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].setInternalStateType(IST_Velocity);
        gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].setIntVarIndx(1);
        gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].setPlotMode(OGC_nodeVectorPlot);
        if ( gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ].getActiveStep() != -1 ) {
            //setupData (gc[OOFEG_DEFORMED_GEOMETRY_LAYER]);
            drawData(gc [ OOFEG_DEFORMED_GEOMETRY_LAYER ]);
        }
    } else if ( mode == pressureField ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_Pressure);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);
        if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
            setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
            drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
        }
    } else if ( mode == vofField ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_VOFFraction);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);
        if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
            setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
            drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
        }
    } else if ( mode == densityField ) {
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_Density);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);
        if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
            setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
            drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
        }
    } else if ( mode == relativeMeshSizeDensity ) {
        //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setInternalStateType(IST_RelMeshDensity);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setIntVarIndx(1);
        gc [ OOFEG_VARPLOT_PATTERN_LAYER ].setPlotMode(OGC_scalarPlot);
        if ( gc [ OOFEG_VARPLOT_PATTERN_LAYER ].getActiveStep() != -1 ) {
            setupData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
            drawData(gc [ OOFEG_VARPLOT_PATTERN_LAYER ]);
        }
    }
}


void
setupSmoother(oofegGraphicContext &gc)
{
    int id;
    FloatArray var;
    IntArray map;

    if ( gc.getActiveStep() == -1 ) {
        return;
    }

    if ( gc.giveIntVarMode() == ISM_recovered ) {
        NodalRecoveryModel *smoother;
        Domain *domain;
        for ( id = 1; id <= gc.getActiveProblem()->giveNumberOfDomains(); id++ ) {
            domain = gc.getActiveProblem()->giveDomain(id);
            smoother = domain->giveSmoother();
            if ( !smoother ) {
                //setup smoother
                setSmoother( gc.giveSmootherType() );
                smoother = domain->giveSmoother();
            }

            smoother->recoverValues( gc.giveIntVarType(), gc.getActiveProblem()->giveCurrentStep() );
            //   problem->giveSmoother()->giveMinMaxVal(&min,&max);
        }
    }

    ColorSetupFringeTableByMinMax( gc.getFringeTable(), ( FPNum ) gc.getScaleMin(), ( FPNum ) gc.getScaleMax() );
    //deleteGraphics(gc);
    //deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
}

void setSmoother(Widget w, XtPointer ptr, XtPointer call_data)
{
    SmootherType mode;
    mode = * ( ( SmootherType * ) ptr );

    setSmoother(mode);
    updateGraphics();
}

void setSmoother(SmootherType mode)
{
    for ( int id = 1; id <= gc [ 0 ].getActiveProblem()->giveNumberOfDomains(); id++ ) {
        if ( mode == Smoother_NA ) {
            gc [ 0 ].getActiveProblem()->giveDomain(id)->
            setSmoother( new NodalAveragingRecoveryModel( gc [ 0 ].getActiveProblem()->giveDomain(id) ) );
        } else if ( mode == Smoother_ZZ ) {
            gc [ 0 ].getActiveProblem()->giveDomain(id)->setSmoother( new ZZNodalRecoveryModel( gc [ 0 ].getActiveProblem()->giveDomain(id) ) );
        } else if ( mode == Smoother_SPR ) {
            gc [ 0 ].getActiveProblem()->giveDomain(id)->setSmoother( new SPRNodalRecoveryModel( gc [ 0 ].getActiveProblem()->giveDomain(id) ) );
        } else {
            OOFEM_ERROR("Unrecognized nodal recovery model");
        }

        gc [ 0 ].getActiveProblem()->giveDomain(id)->giveSmoother()->setRecoveryMode( -1, IntArray() );
    }

    gc [ 0 ].setSmootherType(mode);
}


/*
 * void stressPlot (int i)
 * {
 * DrawMode  mode;
 * if ((i<1)||(i>14)) i = 1;
 * mode  = (DrawMode) (i + internalStateBegin);
 *
 * if (gc.getActiveStep() != -1) {
 *  double min,max;
 *  problem->giveConnectivityTable()->giveMinMaxVal(mode,&min,&max);
 *  if (fabs(max - min) < 1.e-6 ) max += 1.e-6;
 *  ColorSetupFringeTableByMinMax(ft,min,max);
 *  deleteLayerGraphics(OOFEG_STRESS_CONTOUR_LAYER);
 * gc.setDrawMode (mode);
 *  drawData(gc);
 * }
 * }
 */

char *readSimpleString(char *source, char *simpleString, char **remain)
// reads Simple string from source according to following rules:
// at begining skips whitespace (blank, tab)
// read string terminated by whitespace or end-of-line
// remain is unread remain of source string.
{
    char *curr = source;
    char *ss = simpleString;

    while ( ( * curr == ' ' ) || ( * curr == '\n' ) || ( * curr == '\t' ) || ! * curr ) {
        curr++;
    }

    if ( !curr ) {
        fprintf(stderr, "End-of-line encountered\a\n");
        exit(1);
    }

    while ( !( ( * curr == ' ' ) || ( * curr == '\n' ) || ( * curr == '\t' ) || ! * curr ) ) {
        * ss++ = * curr++;
    }

    * ss = '\0';
    * remain = curr;
    return simpleString;
}

/*
 * static LIST toDelete = NULL;
 * int  mark ( NODE data, NODE g_p)
 * {
 * //int aLayer = (int)data;
 * int aLayer = (int)(*data);
 *
 * if (toDelete == NULL)
 *    toDelete = make_list();
 *
 * if (EGGetLayer((GraphicObj *)g_p) == aLayer)
 *    add_to_tail(toDelete, g_p);
 *
 * return 1;
 * }
 *
 *
 * void deleteLayerGraphics (int iLayer)
 * {
 * GraphicObj *g_p;
 *
 * //fprintf (stderr, "(Deleting layer %d ...", iLayer);
 * //EMDispatchToGraphicObjects(age_model, mark, (caddr_t)iLayer);
 * EMDispatchToGraphicObjects(age_model, mark, (caddr_t)&iLayer);
 * if (toDelete == NULL) return ;
 * loop_through_list(toDelete, g_p, GraphicObj *)
 *  EMDeleteGraphics(age_model, g_p);
 * flush_list(toDelete, AND_NODES);
 * //fprintf (stderr, ") done\n");
 * }
 */


void deleteLayerGraphics(int iLayer)
{
    GraphicObj *p;
    int total, deleted;
    EView *v;

    fprintf(stderr, "(Deleting layer %d ...", iLayer);

    v = ( EView * ) get_list_next(age_model->dependent_views, NULL);
    while ( v != NULL ) {
        EVFlushAllDisplayLists(v); // flushes display lists, but keeps nodes
        v = ( EView * ) get_list_next(age_model->dependent_views, v);
    }


    p = ( GraphicObj * ) get_list_next(age_model->graphics_list, NULL);
    total = deleted = 0;
    while ( p != NULL ) {
        total++;
        if ( EGGetLayer(p) == iLayer ) {
            p = ( GraphicObj * ) get_from_current_next(age_model->graphics_list);
            EGDeepDestroyGraphics(p); //free (p);
            p = ( GraphicObj * ) get_list_current(age_model->graphics_list);
            deleted++;
        } else {
            p = ( GraphicObj * ) get_list_next(age_model->graphics_list, p);
        }
    }

    fprintf(stderr, "%d/%d) done\n", total, deleted);

    v = ( EView * ) get_list_next(age_model->dependent_views, NULL);
    while ( v != NULL ) {
        //EMRegenerateGraphics(age_model, v);
        EMRegenerateGraphics(age_model, v); // flushes display lists, but keeps nodes
        //EVFastRedraw(v);
        v = ( EView * ) get_list_next(age_model->dependent_views, v);
    }
}




/************************************************************************/
void print_state(Widget w, XtPointer text_ptr, XtPointer call_data)
/************************************************************************/
{
    EVSetSelectActionFunction(do_print_state, NULL, NULL);
    EMInstallHandler(ESIModel(), EVSelectGraphicsHandler,
                     NULL, NULL,
                     NULL, NULL,
                     NULL, NULL,
                     NULL, NULL);
    ERptStatusMessage(1, 1, "Select element(s)");
}



/****************************************************/
void do_print_state(EView *v_p, caddr_t data_p)
/****************************************************/
/* This is a new version printing to a window. */
{
    GraphicObj *grep = NULL;
    FEMComponent *obj = NULL;
    FILE *fp;
    char buf [ L_tmpnam ];

    if ( EVItemsInSelectionList() <= 0 ) {
        ERptErrMessage(1, 2, "No elements selected.", ERROR_GRADE);
        goto exit;
    }

    tmpnam(buf);
    fp = fopen(buf, "w");
    while ( ( grep = EVGetFromSelectionList() ) != NULL ) {
        if ( EGGetAppObj(grep) != NULL ) {
            // print info
            obj = ( FEMComponent * ) EGGetAppObj(grep);
            if ( gc [ 0 ].getActiveStep() != -1 ) {
                obj->printOutputAt( fp, gc [ 0 ].getActiveProblem()->giveCurrentStep() );
            } else {
                fprintf( fp, "\n%s %d", obj->giveClassName(), obj->giveNumber() );
            }
        }
    }

    fclose(fp);
    ESIPopupInfoFromFile(buf);

exit:
    EVUninstallHandler(v_p);
    ERptStatusMessage(1, 1, "Done");
}


static void apply_layer_update(Widget w, XtPointer ptr, XtPointer call_data)
{
    char tname [ 32 ];
    Widget paned = ( Widget ) ptr;
    Widget t;
    Arg al [ 1 ];
    int ac;
    Boolean s;

    for ( int i = 0; i < OOFEG_LAST_LAYER; i++ ) {
        sprintf(tname, "layer_toggle_%d", i);
        t = XtNameToWidget(paned, tname);
        if ( t != NULL ) {
            ac = 0;
            XtSetArg(al [ ac ], ( String ) XtNstate, & s);
            ac++;
            XtGetValues(t, al, ac);
            gc [ i ].setActivityFlag(s);
        }
    }
}


static void apply_layer_settings(Widget w, XtPointer ptr, XtPointer call_data)
{
    EVSetApplyToViewFunction(set_layer_on_off, ( caddr_t ) ptr, ( ResetFuncP ) uninstall_apply_to_view);
    EMPushHandler(ESIModel(), EVApplyToViewHandler, NULL);
}



static void
set_layer_on_off(EView *v_p, caddr_t data, WCRec *p)
{
    char tname [ 32 ];
    Widget paned = ( Widget ) data;
    Widget t;
    int i;
    Arg al [ 1 ];
    int ac;
    Boolean s;

    int layers [ OOFEG_LAST_LAYER ], dim = 0;
    BOOLEAN on [ OOFEG_LAST_LAYER ];

    for ( i = 0; i < OOFEG_LAST_LAYER; i++ ) {
        sprintf(tname, "layer_toggle_%d", i);
        t = XtNameToWidget(paned, tname);
        if ( t != NULL ) {
            ac = 0;
            XtSetArg(al [ ac ], ( String ) XtNstate, & s);
            ac++;
            XtGetValues(t, al, ac);
            /*      EVSetLayerOnOff(v_p, i, (BOOLEAN)s);  */
            if ( v_p->view_layers [ i ].layer_on != ( BOOLEAN ) s ) {
                layers [ dim ] = i;
                on [ dim ] = ( BOOLEAN ) s;
                dim++;
            }
        }
    }

    if ( dim > 0 ) {
        EVFastSetMultLayersOnOff(v_p, layers, on, dim);
    }
}


static void
uninstall_apply_to_view(EView *v_p, caddr_t data)
{
    EMUninstallHandler( EVGetModel(v_p) );
}

static void apply_mat_reg_filter(Widget data, XtPointer ptr, XtPointer call_data)
{
    char tname [ 32 ];
    Widget paned = ( Widget ) ptr;
    Widget t;
    int i, ii, nmat;
    Arg al [ 1 ];
    int ac;
    Boolean s;


    ii = 0;
    for ( int id = 1; id <= gc [ 0 ].getActiveProblem()->giveNumberOfDomains(); id++ ) {
        nmat = gc [ 0 ].getActiveProblem()->giveDomain(id)->giveNumberOfMaterialModels();

        // int layers[OOFEG_layer_num], dim = 0;
        // BOOLEAN on[OOFEG_layer_num];

        for ( i = 1; i <= nmat; i++ ) {
            ii++;
            sprintf(tname, "mat_reg_toggle_%d-%d", id, i);
            t = XtNameToWidget(paned, tname);
            if ( t != NULL ) {
                ac = 0;
                XtSetArg(al [ ac ], ( String ) XtNstate, & s);
                ac++;
                XtGetValues(t, al, ac);
                /*      EVSetLayerOnOff(v_p, i, (BOOLEAN)s);  */
                gc [ 0 ].setMaterialModelFilterState(ii, ( BOOLEAN ) s);
                //printf ("mat reg filter %d set to %d\n",i,(BOOLEAN)s);
            }
        }
    }
}

static void
pass_setscale_command(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;
    Arg al [ 2 ];
    char *s;
    double max, min;

    ac = 0;
    XtSetArg(al [ ac ], XtNstring, & s);
    ac++;
    XtGetValues(scale_min, al, ac);
    TypeInParseLine(s);
    if ( TypeInGetTokenType(1) != NUMBER ) {
        goto err_exit;
    }

    min = strtod(s, NULL);

    XtSetArg(al [ ac ], XtNstring, & s);
    ac++;
    XtGetValues(scale_max, al, ac);
    TypeInParseLine(s);
    if ( TypeInGetTokenType(1) != NUMBER ) {
        goto err_exit;
    }

    max = strtod(s, NULL);
    if ( max <= min ) {
        printf("setscale_command invalid min. max\n");
        return;
    }


    // gc.setDrawMode (mode);
    // problem->giveConnectivityTable()->giveMinMaxVal(gc,&min,&max);

#if 0
 #define COLOR_SCALE_NUM_LABELS 11
    delta = ( max - min ) / ( COLOR_SCALE_NUM_LABELS - 2 );
    max += delta;
    min -= delta;
#endif

    if ( fabs(max - min) < 1.e-20 ) {
        max += 1.e-20;
    }

    //ColorSetupFringeTableByMinMax(ft,min,max);
    deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
    gc [ 0 ].setScaleVals(min, max);
    gc [ 0 ].setScaleMode(SM_UserDefined);
    updateGraphics();

    return;

err_exit:
    ERptErrMessage(1, 1, "Set scale - bad input; must get 2 numbers", ERROR_GRADE);
}

static void
pass_setanimate_command(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;
    Arg al [ 2 ];
    char *s;
    int estep, sstep;
    int stepinfo [ 2 ];
    int istep, iversion = 0;
    contextIOResultType result;


    ac = 0;
    XtSetArg(al [ ac ], XtNstring, & s);
    ac++;
    XtGetValues(start_step, al, ac);
    TypeInParseLine(s);
    if ( TypeInGetTokenType(1) != NUMBER ) {
        goto err_exit;
    }

    sstep = atoi(s);

    XtSetArg(al [ ac ], XtNstring, & s);
    ac++;
    XtGetValues(end_step, al, ac);
    TypeInParseLine(s);
    if ( TypeInGetTokenType(1) != NUMBER ) {
        goto err_exit;
    }

    estep = atoi(s);
    if ( estep <= sstep ) {
        printf("setanimate_command invalid start/end\n");
        return;
    }

    for ( istep = sstep; istep <= estep; istep++ ) {
        stepinfo [ 0 ] = istep;
        stepinfo [ 1 ] = iversion;
        try {
            result = problem->restoreContext(NULL, CM_State, ( void * ) stepinfo);
        } catch(ContextIOERR & m) {
            m.print();
            return;
        }
        gc [ 0 ].setActiveStep(istep);
        updateGraphics();
    }

    return;

err_exit:
    ERptErrMessage(1, 1, "Set animate - bad input; must get 2 numbers", ERROR_GRADE);
}


static void
pass_setgreyscale_command(Widget w, XtPointer ptr, XtPointer call_data)
{
    int ac;
    Arg al [ 2 ];
    char *s;
    double max, min;

    ac = 0;
    XtSetArg(al [ ac ], XtNstring, & s);
    ac++;
    XtGetValues(greyscale_min, al, ac);
    TypeInParseLine(s);
    if ( TypeInGetTokenType(1) != NUMBER ) {
        goto err_exit;
    }

    min = strtod(s, NULL);

    XtSetArg(al [ ac ], XtNstring, & s);
    ac++;
    XtGetValues(greyscale_max, al, ac);
    TypeInParseLine(s);
    if ( TypeInGetTokenType(1) != NUMBER ) {
        goto err_exit;
    }

    max = strtod(s, NULL);
    if ( max <= min ) {
        printf("setscale_command invalid min. max\n");
        return;
    }


    // gc.setDrawMode (mode);
    // problem->giveConnectivityTable()->giveMinMaxVal(gc,&min,&max);

    if ( fabs(min) < 0.0 ) {
        min = 0.0;
    }

    if ( ( max  < 0.0 ) || ( max > 1.0 ) ) {
        max = 1.0;
    }

    GreyValueMinMaxSetup(min, max);
    updateGraphics();

    return;

err_exit:
    ERptErrMessage(1, 1, "Set scale - bad input; must get 2 numbers", ERROR_GRADE);
}

static void
OOFEGReturnHitInCmd(Widget w, XEvent *event, String *params,
                    Cardinal *num_params)
{
    /*
     * Arg args[1];
     * String str;
     *
     * XtSetArg(args[0], XtNstring, &str);
     * XtGetValues(w, args, 1);
     *
     * // exec_cmd(str);
     * printf ("%s\n",str);
     *
     * XtSetArg(args[0], XtNstring, str);
     * XtSetValues(w, args, 1);
     */
}


static void apply_toggleIcs(Widget w, XtPointer ptr, XtPointer call_data)
{
    EVSetApplyToViewFunction(toggleIcs, ( caddr_t ) ptr, ( ResetFuncP ) uninstall_apply_to_view);
    EMPushHandler(ESIModel(), EVApplyToViewHandler, NULL);
}

static void
toggleIcs(EView *v_p, caddr_t data, WCRec *p)
{
    EVToggleScale(v_p);
}

static void
toggle_IcsColors(Widget w, XtPointer ptr, XtPointer call_data)
{
    RevertFringeTableColorScale();
}

static void
setIcsToColor(Widget w, XtPointer ptr, XtPointer call_data)
{
    SwitchFringeTableToColor();
}

static void
setIcsToGrey(Widget w, XtPointer ptr, XtPointer call_data)
{
    SwitchFringeTableToGrey();
}

static void
toggleSmoothScale(Widget w, XtPointer ptr, XtPointer call_data)
{
    ToggleFringeTableSmoothScale();
}

static void
toggleTransparentContours(Widget w, XtPointer ptr, XtPointer call_data)
{
    if ( IsContourBgTransparent() ) {
        SetContourBgTransparent(true);
    } else {
        SetContourBgTransparent(true);
    }
}





static void
setAutoscale(Widget w, XtPointer ptr, XtPointer call_data)
{
    gc [ 0 ].setScaleMode(SM_Autoscale);
    updateGraphics();
}

void
deleteGraphics(oofegGraphicContext &gc)
{
    OGC_PlotModeType type = gc.giveIntVarPlotMode();
    switch ( type ) {
    case OGC_rawGeometry:
        deleteLayerGraphics(OOFEG_RAW_GEOMETRY_LAYER);
        break;
    case OGC_deformedGeometry:
    case OGC_eigenVectorGeometry:
    case OGC_nodeVectorPlot:
        deleteLayerGraphics(OOFEG_DEFORMED_GEOMETRY_LAYER);
        break;
    case OGC_nodeAnnotation:
    case OGC_nodeGeometry:
        deleteLayerGraphics(OOFEG_NODE_ANNOTATION_LAYER);
        break;
    case OGC_elementAnnotation:
        deleteLayerGraphics(OOFEG_ELEMENT_ANNOTATION_LAYER);
        break;
    case OGC_essentialBC:
        deleteLayerGraphics(OOFEG_BCIC_ANNOTATION_LAYER);
        break;
    case OGC_naturalBC:
        deleteLayerGraphics(OOFEG_NATURALBC_LAYER);
        break;
    case OGC_scalarPlot:
    case OGC_vectorPlot:
    case OGC_tensorPlot:
        deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
        break;
    case OGC_elemSpecial:
        deleteLayerGraphics(OOFEG_CRACK_PATTERN_LAYER);
        break;
    default:
        deleteLayerGraphics(OOFEG_VARPLOT_PATTERN_LAYER);
    }
}


void setupData(oofegGraphicContext &gc)
{
    gc.resetScaleVals();
    setupSmoother(gc);

    // if isoline plot mode is active then determine the scale min and max values from the smoother
    // because the color scale must be properly set up before isolines are drawn, since the
    // isolines are drawn at positions corresponding to values of color scale.
    // Note: isolines not allowed on local values.
    /*
     * if (gc.getScalarAlgo() == SA_ISO_LINE) {
     * int i,j,enodes,esides,nelems, indx;
     * int init = 1;
     * double minVal, maxVal, stress;
     * Element *elem;
     * FloatArray var;
     * IntArray map;
     *
     * Domain* domain = gc.getActiveProblem()->giveDomain(gc.getActiveDomain());
     * nelems = domain-> giveNumberOfElements();
     * for (i = 1; i<= nelems; i++) {
     * elem = domain -> giveElement(i);
     * enodes = elem ->giveNumberOfNodes();
     * // loop over each element node
     * for (j = 1; j<= enodes; j++) {
     * if (elem -> giveInternalStateAtNode (var, gc.giveIntVarType(), gc.giveIntVarMode(),j,problem->giveCurrentStep())) {
     *  elem->giveIntVarCompFullIndx (map, gc.giveIntVarType());
     *  if ((indx = map.at(gc.giveIntVarIndx()))) {
     *   stress = var.at(indx);
     *   if (init) {maxVal = minVal = stress; init = 0;}
     *   if (stress > maxVal) maxVal = stress;
     *   if (stress < minVal) minVal = stress;
     *  }
     * }
     * }
     *
     * esides = elem ->giveNumberOfSides();
     * // loop over each element side
     * for (j = 1; j<= esides; j++) {
     *  if (elem -> giveInternalStateAtSide (var, gc.giveIntVarType(), gc.giveIntVarMode(),j,problem->giveCurrentStep())) {
     *   elem->giveIntVarCompFullIndx (map, gc.giveIntVarType());
     *   if ((indx = map.at(gc.giveIntVarIndx()))) {
     *    stress = var.at(indx);
     *    if (init) {maxVal = minVal = stress; init = 0;}
     *    if (stress > maxVal) maxVal = stress;
     *    if (stress < minVal) minVal = stress;
     *   }
     *  }
     * }
     * }
     *
     * ColorSetupFringeTableByMinMax(gc.getFringeTable(),(FPNum) minVal,(FPNum) maxVal);
     * }
     */
}


void showSparseMtrxStructure(Widget wid, XtPointer cl, XtPointer cd)
{
    if ( gc [ OOFEG_SPARSE_PROFILE_LAYER ].getActiveStep() != -1 ) {
        TimeStep *tStep = gc [ OOFEG_SPARSE_PROFILE_LAYER ].getActiveProblem()->giveCurrentStep();
        gc [ OOFEG_SPARSE_PROFILE_LAYER ].setSparseProfileMode(1);
        deleteLayerGraphics(OOFEG_SPARSE_PROFILE_LAYER);
        gc [ OOFEG_SPARSE_PROFILE_LAYER ].getActiveProblem()->showSparseMtrxStructure(1, gc [ OOFEG_SPARSE_PROFILE_LAYER ], tStep);
    }
}




static Widget
oofeg_add_palette(const char *palette_label, Widget parent, Widget *palette)
{
    Widget button;

    /* create palette with simple palette label and then apply modifications */

    button = ESIAddPalette(const_cast< char * >(OOFEG_PALETTE_BUTTON_RESOURCE), const_cast< char * >(palette_label), parent, NULL, 0, palette);
    return ( button );
}

static Widget
oofeg_add_popdown_menu(const char *menu_label, Widget parent, Widget *palette)
{
    return ( ESIAddPopdownMenu(const_cast< char * >(OOFEG_MENU_BUTTON_RESOURCE), const_cast< char * >(menu_label), parent, NULL, 0, palette) );
}

static Widget
oofeg_add_button(const char *name, const char *button_label, WidgetClass wclass, Widget palette,
                 XtCallbackProc action, XtPointer data)
{
    return ( ESIAddButton(const_cast< char * >(name), const_cast< char * >(button_label), wclass, palette, NULL, 0, action, data) );
}

static Widget
oofeg_add_button(const char *name, const char *button_label, WidgetClass wclass, Widget palette, Arg *arg, int ac,
                 XtCallbackProc action, XtPointer data)
{
    return ( ESIAddButton(const_cast< char * >(name), const_cast< char * >(button_label), wclass, palette, arg, ac, action, data) );
}

static Widget
oofeg_add_menu_item(const char *name, const char *item_label, Widget palette,
                    XtCallbackProc action, XtPointer data)
{
    return ( ESIAddMenuItem(const_cast< char * >(name), const_cast< char * >(item_label), palette, NULL, 0, action, data) );
}



static Widget
oofeg_add_dialog(const char *name, const char *dialog_label, const char *prompt, const char *init_value, Widget palette,
                 XtCallbackProc action, const char *data, ESIDialogValueType type, ESIVerifyValueProc proc)
{
    Widget button;
    button = ESIAddPopupDialog(const_cast< char * >(name), const_cast< char * >(dialog_label), const_cast< char * >(prompt),
                               const_cast< char * >(init_value), palette,
                               NULL, 0, action, ( XtPointer ) data, type, proc);

    return ( button );
}


void
oofeg_display_message(const char *message)
{
    /* setting message to NULL causes display of last previous message */

    ERptStatusMessage(1, 1, message);
}

void
oofeg_exit(Widget w, XtPointer ptr, XtPointer call_data)
{
    ESIPopupConfirmDialog( ( Widget ) ptr, const_cast< char * >("Really exit?"), oofeg_quit, ptr, NULL, NULL );
}

void
oofeg_quit(XtPointer ptr)
{
    XtCloseDisplay( XtDisplay( ESITopLevelWidget() ) );
    exit( ESISetExitCode(0) );
}

void
oofeg_open_frame(Widget w, XtPointer ptr, XtPointer call_data)
{
    EView *view;
    view = ElixirNewView(viewTitle.c_str(), const_cast< char * >("SimpleXF"), const_cast< char * >(OOFEG_BACKGROUND_COLOR),
                         const_cast< char * >(OOFEG_DEFAULTDRAW_COLOR), 500, 400);

    EMAttachView(ESIModel(), view);
    EMRegenerateGraphics(ESIModel(), view);
    EVShowAxes(view, oofeg_axes);
}


void
oofeg_close_frame(Widget w, XtPointer ptr, XtPointer call_data)
{
    EVSetApplyToViewFunction( ( void(*) (EView *, caddr_t, WCRec *) )oofeg_destroy_frame, NULL, NULL );
    EVSetApplyToViewPreventRedisplay(YES);
    EMPushHandler(ESIModel(), EVApplyToViewHandler, NULL);
}

void
oofeg_destroy_frame(EView *view, caddr_t data, WCRec *p)
{
    EVFreeView(view);
    XtDestroyWidget(view->top_view_widget);
}

void
oofeg_set_render(Widget w, XtPointer ptr, XtPointer call_data)
{
    EVSetApplyToViewFunction(set_render, ( caddr_t ) ptr, ( ResetFuncP ) uninstall_apply_to_view);
    EMPushHandler(ESIModel(), EVApplyToViewHandler, NULL);
}

void
set_render(EView *v_p, caddr_t data, WCRec *p) {
    EView *view = v_p;
    int rt = * ( ( int * ) data );

    switch ( rt ) {
    case OOFEG_WIRED_RENDER:
        if ( EVGetRenderMode(view) != WIRE_RENDERING ) {
            EVSetRenderMode(view, WIRE_RENDERING);
            EVFastRedraw(view);
        }

        break;
    case OOFEG_NORMAL_RENDER:
        if ( EVGetRenderMode(view) != NORMAL_RENDERING ) {
            EVSetRenderMode(view, NORMAL_RENDERING);
            EVFastRedraw(view);
        }

        break;
    case OOFEG_HIDDEN_RENDER:
        if ( EVGetRenderMode(view) != FILLED_HIDDEN_RENDERING ) {
            EVSetRenderMode(view, FILLED_HIDDEN_RENDERING);
            EVFastRedraw(view);
        }

        break;
    case OOFEG_SHADED_RENDER:
        if ( EVGetRenderMode(view) != CONST_SHADING_RENDERING || EVGetShadeMode(view) != COLOR_SHADING ) {
            EVSetShadeMode(view, COLOR_SHADING);
            EVSetRenderMode(view, CONST_SHADING_RENDERING);
            EVFastRedraw(view);
        }

        break;
    default:
        break;
    }
}

void oofeg_toggle_axes(Widget w, XtPointer ptr, XtPointer call_data)
{
    if ( oofeg_axes == 0 ) {
        EMDispatchToDependentViews(ESIModel(), show_axes, NULL);
        oofeg_axes = 1;
    } else {
        EMDispatchToDependentViews(ESIModel(), hide_axes, NULL);
        oofeg_axes = 0;
    }
}



int
show_axes(NODE _data, NODE _view)
{
    EVShowAxes( ( EView * ) _view, YES );
    return ( 1 );
}

int
hide_axes(NODE _data, NODE _view)
{
    EVShowAxes( ( EView * ) _view, NO );
    return ( 1 );
}


void drawData(oofegGraphicContext &gc)
{
    // switch some off render modes
    EView *v;
    IntArray rm( list_length ( age_model->dependent_views ) );
    ERenderingType rt;
    int i = 1;

    if ( gc.giveIntVarPlotMode() == OGC_unknown ) {
        return;
    }

    v = ( EView * ) get_list_next(age_model->dependent_views, NULL);
    while ( v != NULL ) {
        rt = EVGetRenderMode(v);
        rm.at(i++) = ( int ) rt;
        if ( ( rt == CONST_SHADING_RENDERING ) || ( rt == FILLED_HIDDEN_RENDERING ) ) {
            EVSetRenderMode(v, NORMAL_RENDERING);
        }

        v = ( EView * ) get_list_next(age_model->dependent_views, v);
    }

    deleteGraphics(gc);
    gc.getActiveProblem()->drawYourself(gc);

    // restore orinal rendering modes
    i = 1;
    v = ( EView * ) get_list_next(age_model->dependent_views, NULL);
    while ( v != NULL ) {
        rt = EVGetRenderMode(v);
        if ( ( ( int ) rt ) != rm.at(i) ) {
            EVSetRenderMode( v, ( ERenderingType ) rm.at(i) );
        }

        EVFastRedraw(v);
        i++;
        v = ( EView * ) get_list_next(age_model->dependent_views, v);
    }

    //printf ("Color Scale (%e, %e)\n", gc.getScaleMin(),gc.getScaleMax());
}

int
fit_all(NODE _data, NODE _view)
{
    EVFitAllIntoView( ( EView * ) _view );
    return ( 1 );
}


void
oofeg_fit_all_graphics(Widget w, XtPointer ptr, XtPointer call_data)
{
    EMDispatchToDependentViews(ESIModel(), fit_all, NULL);
}

static void
set_bg_color(Widget w, XtPointer ptr, XtPointer call_data)
{
    BOOLEAN success;
    char *bg_color = ( char * ) ptr, oofeg_bg_color [ 64 ];

    strncpy(oofeg_bg_color, bg_color, 64);
    oofeg_background_color = ColorGetPixelFromString(oofeg_bg_color, & success);
    if ( success == NO ) {
        fprintf(stderr, "Unable to allocate %s as background color", oofeg_bg_color);
        return;
    }

    EMDispatchToDependentViews(ESIModel(), set_background, NULL);
}

static int
set_background(NODE _data, NODE _view)
{
    EVSetBackground( ( EView * ) _view, oofeg_background_color );
    return ( 1 );
}

void
updateISA(oofegGraphicContext *context)
{
    char buff [ 1024 ];
    double time;
    TimeStep *ct = context [ 0 ].getActiveProblem()->giveCurrentStep();
    if ( ct != NULL ) {
        time = ct->giveTargetTime();
    } else {
        time = 0.0;
    }

    if ( GetISASize() > 2 ) {
        sprintf(buff, "| %s | step: %4d.%d | time: %e |", jobName.c_str(), context [ 0 ].getActiveStep(), context [ 0 ].getActiveStepVersion(),
                time);
        SetISAContent(1, buff);
    }
}


static void toggleStatus(Widget w, XtPointer ptr, XtPointer call_data)
{
    EVSetApplyToViewFunction(view_toggleStatus_on_off, ( caddr_t ) ptr, ( ResetFuncP ) uninstall_apply_to_view);
    EMPushHandler(ESIModel(), EVApplyToViewHandler, NULL);
}

static void
view_toggleStatus_on_off(EView *v_p, caddr_t data, WCRec *p)
{
    Boolean status = v_p->show_status;
    if ( status == NO ) {
        EVShowStatus(v_p, YES);
    } else {
        EVShowStatus(v_p, NO);
    }
}

void updateGraphics()
{
    updateISA(gc);
    for ( int i = 0; i < OOFEG_LAST_LAYER; i++ ) {
        if ( gc [ i ].isActive() ) {
            if ( gc [ i ].giveIntVarPlotMode() != OGC_unknown ) {
                if ( i == OOFEG_VARPLOT_PATTERN_LAYER ) {
                    setupData(gc [ i ]);                       // setup data only for varplots
                }

                drawData(gc [ i ]);
            }

            //printf ("Color Scale_%d (%e, %e)\n", i, gc[i].getScaleMin(),gc[i].getScaleMax());
        }
    }
}

#ifdef OOFEG_DEVEL
void debug_run(Widget w, XtPointer ptr, XtPointer call_data)
{
    try {
        problem->checkProblemConsistency();
        gc [ 0 ].setActiveStep(0);
        problem->solveYourself();
    } catch(OOFEM_Terminate & c) {
        delete problem;
 #ifdef __PETSC_MODULE
        PetscFinalize();
 #endif
 #ifdef __PARALLEL_MODE
        MPI_Finalize();
 #endif
    }
}
#endif
#endif // __OOFEG
