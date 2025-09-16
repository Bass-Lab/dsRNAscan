/* Minimal EMBOSS compatibility header for standalone einverted
 * Implements just enough EMBOSS functions to compile einverted.c
 */

#ifndef EMBOSS_MINIMAL_H
#define EMBOSS_MINIMAL_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdarg.h>
#include <limits.h>

/* Basic types */
typedef int ajint;
typedef unsigned int ajuint;
typedef long ajlong;
typedef int AjBool;
typedef enum {ajEOutfileReturnUnknown} AjEOutfileReturn;

#define ajTrue 1
#define ajFalse 0
#define AJMAX(a,b) ((a)>(b)?(a):(b))
#define AJMIN(a,b) ((a)<(b)?(a):(b))

/* String structure */
typedef struct {
    char *Ptr;
    size_t Len;
    size_t Res;
    ajuint Use;
} AjOStr;
typedef AjOStr* AjPStr;

/* File structure */
typedef struct {
    FILE *fp;
    AjPStr Name;
} AjSFile;
typedef AjSFile* AjPFile;

/* Sequence structure */
typedef struct {
    AjPStr Name;
    AjPStr Seq;
    ajint Begin;
    ajint End;
    ajint Length;
    ajint Offset;
    ajint Offend;
    AjBool Rev;
    AjBool Reversed;
} AjOSeq;
typedef AjOSeq* AjPSeq;

/* Sequence input structure */
typedef struct {
    FILE *File;
    AjPSeq Seq;
} AjOSeqall;
typedef AjOSeqall* AjPSeqall;

/* Sequence output structure */
typedef struct {
    AjPFile File;
    AjPStr Name;
} AjOSeqout;
typedef AjOSeqout* AjPSeqout;

/* Integer array structure */
typedef struct {
    ajint *Ptr;
    ajuint Len;
    ajuint Res;
} AjOInt;
typedef AjOInt* AjPInt;

/* 2D integer array structure */
typedef struct {
    ajint **Ptr;
    ajuint Len;
    ajuint Res;
} AjOInt2d;
typedef AjOInt2d* AjPInt2d;

/* List structure (simplified) */
typedef struct {
    void **items;
    ajuint count;
    ajuint size;
} AjOList;
typedef AjOList* AjPList;

/* Sequence converter */
typedef struct {
    char *chars;
} AjOSeqCvt;
typedef AjOSeqCvt* AjPSeqCvt;

/* Debug macro - disable for production */
#define ajDebug if(0) printf
#define ajUser printf

/* Function implementations */
static inline AjPStr ajStrNew(void) {
    AjPStr str = (AjPStr)calloc(1, sizeof(AjOStr));
    str->Ptr = (char*)malloc(1);
    str->Ptr[0] = '\0';
    str->Len = 0;
    str->Res = 1;
    str->Use = 1;
    return str;
}

static inline AjPStr ajStrNewC(const char *txt) {
    AjPStr str = (AjPStr)calloc(1, sizeof(AjOStr));
    str->Len = strlen(txt);
    str->Res = str->Len + 1;
    str->Ptr = (char*)malloc(str->Res);
    strcpy(str->Ptr, txt);
    str->Use = 1;
    return str;
}

static inline void ajStrDel(AjPStr *pstr) {
    if (!pstr || !*pstr) return;
    free((*pstr)->Ptr);
    free(*pstr);
    *pstr = NULL;
}

static inline void ajStrAssignC(AjPStr *pstr, const char *txt) {
    if (!*pstr) *pstr = ajStrNew();
    size_t len = strlen(txt);
    if ((*pstr)->Res <= len) {
        (*pstr)->Res = len + 1;
        (*pstr)->Ptr = (char*)realloc((*pstr)->Ptr, (*pstr)->Res);
    }
    strcpy((*pstr)->Ptr, txt);
    (*pstr)->Len = len;
}

static inline const char* ajStrGetPtr(const AjPStr str) {
    return str ? str->Ptr : "";
}

static inline ajuint ajStrGetLen(const AjPStr str) {
    return str ? str->Len : 0;
}

static inline AjPFile ajFileNewOutNameC(const char *name) {
    AjPFile file = (AjPFile)calloc(1, sizeof(AjSFile));
    if (strcmp(name, "stdout") == 0) {
        file->fp = stdout;
    } else {
        file->fp = fopen(name, "w");
    }
    file->Name = ajStrNewC(name);
    return file;
}

static inline void ajFileClose(AjPFile *pfile) {
    if (!pfile || !*pfile) return;
    if ((*pfile)->fp && (*pfile)->fp != stdout && (*pfile)->fp != stderr) {
        fclose((*pfile)->fp);
    }
    ajStrDel(&(*pfile)->Name);
    free(*pfile);
    *pfile = NULL;
}

static inline void ajFmtPrintF(AjPFile file, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    vfprintf(file ? file->fp : stdout, fmt, args);
    va_end(args);
}

static inline AjPSeq ajSeqNew(void) {
    AjPSeq seq = (AjPSeq)calloc(1, sizeof(AjOSeq));
    seq->Name = ajStrNew();
    seq->Seq = ajStrNew();
    return seq;
}

static inline void ajSeqDel(AjPSeq *pseq) {
    if (!pseq || !*pseq) return;
    ajStrDel(&(*pseq)->Name);
    ajStrDel(&(*pseq)->Seq);
    free(*pseq);
    *pseq = NULL;
}

static inline ajint ajSeqGetLen(const AjPSeq seq) {
    return seq ? seq->Length : 0;
}

static inline const AjPStr ajSeqGetNameS(const AjPSeq seq) {
    return seq ? seq->Name : NULL;
}

static inline void ajSeqAssignNameS(AjPSeq seq, const AjPStr name) {
    if (seq && name) {
        if (!seq->Name) seq->Name = ajStrNew();
        ajStrAssignC(&seq->Name, ajStrGetPtr(name));
    }
}

static inline void ajSeqAssignSeqS(AjPSeq seq, const AjPStr seqstr) {
    if (seq && seqstr) {
        if (!seq->Seq) seq->Seq = ajStrNew();
        ajStrAssignC(&seq->Seq, ajStrGetPtr(seqstr));
        seq->Length = ajStrGetLen(seqstr);
    }
}

static inline AjBool ajSeqIsReversed(const AjPSeq seq) {
    return seq ? seq->Reversed : ajFalse;
}

static inline ajint ajSeqGetOffset(const AjPSeq seq) {
    return seq ? seq->Offset : 0;
}

static inline ajint ajSeqGetOffend(const AjPSeq seq) {
    return seq ? seq->Offend : 0;
}

static inline void ajSeqTrim(AjPSeq seq) {
    /* Simplified - just ensure sequence is clean */
    if (!seq || !seq->Seq) return;
}

/* Integer array functions */
static inline AjPInt ajIntNew(void) {
    AjPInt arr = (AjPInt)calloc(1, sizeof(AjOInt));
    return arr;
}

static inline AjPInt ajIntNewRes(ajuint size) {
    AjPInt arr = (AjPInt)calloc(1, sizeof(AjOInt));
    arr->Ptr = (ajint*)calloc(size, sizeof(ajint));
    arr->Res = size;
    arr->Len = 0;
    return arr;
}

static inline void ajIntDel(AjPInt *parr) {
    if (!parr || !*parr) return;
    free((*parr)->Ptr);
    free(*parr);
    *parr = NULL;
}

static inline ajint ajIntGet(const AjPInt arr, ajuint pos) {
    if (!arr || pos >= arr->Res) return 0;
    return arr->Ptr[pos];
}

static inline AjBool ajIntPut(AjPInt *parr, ajuint pos, ajint val) {
    if (!parr || !*parr || pos >= (*parr)->Res) return ajFalse;
    (*parr)->Ptr[pos] = val;
    if (pos >= (*parr)->Len) (*parr)->Len = pos + 1;
    return ajTrue;
}

/* 2D integer array functions */
static inline AjPInt2d ajInt2dNewResRes2(ajuint rows, ajuint cols) {
    AjPInt2d arr = (AjPInt2d)calloc(1, sizeof(AjOInt2d));
    arr->Ptr = (ajint**)calloc(rows, sizeof(ajint*));
    for (ajuint i = 0; i < rows; i++) {
        arr->Ptr[i] = (ajint*)calloc(cols, sizeof(ajint));
    }
    arr->Len = rows;
    arr->Res = cols;
    return arr;
}

static inline void ajInt2dDel(AjPInt2d *parr) {
    if (!parr || !*parr) return;
    for (ajuint i = 0; i < (*parr)->Len; i++) {
        free((*parr)->Ptr[i]);
    }
    free((*parr)->Ptr);
    free(*parr);
    *parr = NULL;
}

static inline AjBool ajInt2dPut(AjPInt2d *parr, ajuint row, ajuint col, ajint val) {
    if (!parr || !*parr || row >= (*parr)->Len || col >= (*parr)->Res) return ajFalse;
    (*parr)->Ptr[row][col] = val;
    return ajTrue;
}

/* List functions (simplified) */
static inline AjPList ajListstrNew(void) {
    AjPList list = (AjPList)calloc(1, sizeof(AjOList));
    list->size = 10;
    list->items = (void**)calloc(list->size, sizeof(void*));
    list->count = 0;
    return list;
}

static inline void ajListstrPushAppend(AjPList list, AjPStr str) {
    if (!list) return;
    if (list->count >= list->size) {
        list->size *= 2;
        list->items = (void**)realloc(list->items, list->size * sizeof(void*));
    }
    list->items[list->count++] = str;
}

static inline void ajListstrFreeData(AjPList *plist) {
    if (!plist || !*plist) return;
    for (ajuint i = 0; i < (*plist)->count; i++) {
        ajStrDel((AjPStr*)&(*plist)->items[i]);
    }
    free((*plist)->items);
    free(*plist);
    *plist = NULL;
}

/* Base comparison including G-U wobble */
static inline AjBool ajBaseAlphacharCompare(char base1, char base2) {
    base1 = toupper(base1);
    base2 = toupper(base2);
    
    /* Standard Watson-Crick */
    if ((base1 == 'A' && (base2 == 'T' || base2 == 'U')) ||
        ((base1 == 'T' || base1 == 'U') && base2 == 'A') ||
        (base1 == 'C' && base2 == 'G') ||
        (base1 == 'G' && base2 == 'C')) {
        return ajTrue;
    }
    
    /* G-U wobble pairs */
    if ((base1 == 'G' && (base2 == 'T' || base2 == 'U')) ||
        ((base1 == 'T' || base1 == 'U') && base2 == 'G')) {
        return ajTrue;
    }
    
    return ajFalse;
}

/* Sequence converter */
static inline AjPSeqCvt ajSeqcvtNewEndC(const char *chars) {
    AjPSeqCvt cvt = (AjPSeqCvt)calloc(1, sizeof(AjOSeqCvt));
    cvt->chars = strdup(chars);
    return cvt;
}

static inline void ajSeqcvtDel(AjPSeqCvt *pcvt) {
    if (!pcvt || !*pcvt) return;
    free((*pcvt)->chars);
    free(*pcvt);
    *pcvt = NULL;
}

static inline void ajSeqConvertNum(const AjPSeq seq, const AjPSeqCvt cvt, AjPStr *pstr) {
    if (!seq || !seq->Seq || !pstr) return;
    
    if (!*pstr) *pstr = ajStrNew();
    
    const char *seqstr = ajStrGetPtr(seq->Seq);
    ajint len = ajSeqGetLen(seq);
    
    if ((*pstr)->Res < (size_t)len + 1) {
        (*pstr)->Res = len + 1;
        (*pstr)->Ptr = (char*)realloc((*pstr)->Ptr, (*pstr)->Res);
    }
    
    /* Convert to numeric: A=0, C=1, G=2, T=3, other=4 */
    for (ajint i = 0; i < len; i++) {
        switch(toupper(seqstr[i])) {
            case 'A': (*pstr)->Ptr[i] = 0; break;
            case 'C': (*pstr)->Ptr[i] = 1; break;
            case 'G': (*pstr)->Ptr[i] = 2; break;
            case 'T': case 'U': (*pstr)->Ptr[i] = 3; break;
            default: (*pstr)->Ptr[i] = 4; break;
        }
    }
    (*pstr)->Len = len;
}

/* Missing macros */
#define AJCNEW(p,c) ((p) = calloc((c), sizeof(*p)))
#define AJFREE(p) do { free(p); (p) = NULL; } while(0)

/* Missing functions */
static inline void ajFmtPrintS(AjPStr *pstr, const char *fmt, ...) {
    va_list args;
    va_start(args, fmt);
    char buffer[10000];
    vsnprintf(buffer, sizeof(buffer), fmt, args);
    va_end(args);
    if (!*pstr) *pstr = ajStrNew();
    ajStrAssignC(pstr, buffer);
}

static inline void ajStrAppendK(AjPStr *pstr, char c) {
    if (!*pstr) *pstr = ajStrNew();
    size_t len = (*pstr)->Len;
    if (len + 2 > (*pstr)->Res) {
        (*pstr)->Res = (len + 2) * 2;
        (*pstr)->Ptr = realloc((*pstr)->Ptr, (*pstr)->Res);
    }
    (*pstr)->Ptr[len] = c;
    (*pstr)->Ptr[len + 1] = '\0';
    (*pstr)->Len = len + 1;
}

static inline void ajStrReverse(AjPStr *pstr) {
    if (!*pstr || (*pstr)->Len == 0) return;
    char *start = (*pstr)->Ptr;
    char *end = start + (*pstr)->Len - 1;
    while (start < end) {
        char tmp = *start;
        *start = *end;
        *end = tmp;
        start++;
        end--;
    }
}

static inline AjBool ajListstrPop(AjPList list, AjPStr *pstr) {
    if (!list || list->count == 0) return ajFalse;
    *pstr = (AjPStr)list->items[0];
    for (ajuint i = 1; i < list->count; i++) {
        list->items[i-1] = list->items[i];
    }
    list->count--;
    return ajTrue;
}

static inline void ajListstrFree(AjPList *plist) {
    if (!plist || !*plist) return;
    free((*plist)->items);
    free(*plist);
    *plist = NULL;
}

static inline void ajSeqoutClose(AjPSeqout seqout) {
    /* Do nothing */
}

static inline void embExit(void) {
    /* Do nothing */
}

/* Global variables for command-line arguments */
static int g_argc = 0;
static char **g_argv = NULL;
static FILE *g_input = NULL;
static FILE *g_output = NULL;
static FILE *g_outseq = NULL;
static int g_sbegin = 1;
static int g_send = 0;
static int g_filter = 0;

/* Parse command-line arguments and store them */
static inline void embInit(const char *pgm, int argc, char **argv) {
    g_argc = argc;
    g_argv = argv;
    g_input = stdin;
    g_output = stdout;
    
    /* Parse arguments */
    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "-sequence") == 0 && i + 1 < argc) {
            i++;
            if (strcmp(argv[i], "stdin") != 0) {
                g_input = fopen(argv[i], "r");
                if (!g_input) g_input = stdin;
            }
        } else if (strcmp(argv[i], "-outfile") == 0 && i + 1 < argc) {
            i++;
            if (strcmp(argv[i], "stdout") != 0) {
                g_output = fopen(argv[i], "w");
                if (!g_output) g_output = stdout;
            }
        } else if (strcmp(argv[i], "-outseq") == 0 && i + 1 < argc) {
            i++;
            if (strcmp(argv[i], "/dev/null") == 0 || strcmp(argv[i], "null") == 0) {
                g_outseq = NULL;
            } else {
                g_outseq = fopen(argv[i], "w");
            }
        } else if (strcmp(argv[i], "-sbegin") == 0 && i + 1 < argc) {
            g_sbegin = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-send") == 0 && i + 1 < argc) {
            g_send = atoi(argv[++i]);
        } else if (strcmp(argv[i], "-filter") == 0) {
            g_filter = 1;
        }
    }
}

static inline AjPFile ajAcdGetOutfile(const char *name) {
    AjPFile file = (AjPFile)calloc(1, sizeof(AjSFile));
    file->fp = g_output;
    return file;
}

static inline ajint ajAcdGetInt(const char *name) {
    /* Look for command-line argument */
    for (int i = 1; i < g_argc - 1; i++) {
        char arg[100];
        sprintf(arg, "-%s", name);
        if (strcmp(g_argv[i], arg) == 0) {
            return atoi(g_argv[i + 1]);
        }
    }
    
    /* Return defaults */
    if (strcmp(name, "threshold") == 0) return 50;
    if (strcmp(name, "match") == 0) return 3;
    if (strcmp(name, "mismatch") == 0) return -4;
    if (strcmp(name, "gap") == 0) return 12;
    if (strcmp(name, "maxrepeat") == 0) return 2000;
    return 0;
}

static inline AjPSeqall ajAcdGetSeqall(const char *name) {
    AjPSeqall seqall = (AjPSeqall)calloc(1, sizeof(AjOSeqall));
    seqall->File = g_input;
    seqall->Seq = NULL;
    return seqall;
}

static inline AjPSeqout ajAcdGetSeqout(const char *name) {
    return NULL; /* We won't use sequence output */
}

static inline AjBool ajSeqallNext(AjPSeqall seqall, AjPSeq *pseq) {
    if (!seqall || !pseq) return ajFalse;
    
    char line[10000];
    static char seq_buffer[10000000];
    static char name_buffer[1000];
    static int seq_len = 0;
    static int reading = 0;
    
    while (fgets(line, sizeof(line), seqall->File)) {
        if (line[0] == '>') {
            if (seq_len > 0) {
                /* Return the previous sequence */
                if (!*pseq) *pseq = ajSeqNew();
                ajStrAssignC(&(*pseq)->Name, name_buffer);
                ajStrAssignC(&(*pseq)->Seq, seq_buffer);
                (*pseq)->Length = seq_len;
                
                /* Start new sequence */
                strncpy(name_buffer, line + 1, 999);
                name_buffer[strcspn(name_buffer, "\n\r")] = '\0';
                seq_len = 0;
                reading = 1;
                
                return ajTrue;
            }
            
            /* Start first sequence */
            strncpy(name_buffer, line + 1, 999);
            name_buffer[strcspn(name_buffer, "\n\r")] = '\0';
            seq_len = 0;
            reading = 1;
        } else if (reading) {
            /* Add to sequence */
            int len = strlen(line);
            if (line[len-1] == '\n') line[--len] = '\0';
            for (int i = 0; i < len; i++) {
                if (isalpha(line[i])) {
                    seq_buffer[seq_len++] = toupper(line[i]);
                }
            }
            seq_buffer[seq_len] = '\0';
        }
    }
    
    /* Return last sequence if any */
    if (seq_len > 0) {
        if (!*pseq) *pseq = ajSeqNew();
        ajStrAssignC(&(*pseq)->Name, name_buffer);
        ajStrAssignC(&(*pseq)->Seq, seq_buffer);
        (*pseq)->Length = seq_len;
        seq_len = 0;
        return ajTrue;
    }
    
    return ajFalse;
}

static inline void ajSeqallDel(AjPSeqall *pseqall) {
    if (!pseqall || !*pseqall) return;
    free(*pseqall);
    *pseqall = NULL;
}

static inline void ajSeqoutWriteSeq(AjPSeqout seqout, const AjPSeq seq) {
    /* Do nothing - we're not writing sequences */
}

static inline void ajSeqoutDel(AjPSeqout *pseqout) {
    /* Do nothing */
}

#endif /* EMBOSS_MINIMAL_H */