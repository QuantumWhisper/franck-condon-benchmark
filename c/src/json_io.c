#include "json_io.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static char *read_file_text(const char *filepath) {
    FILE *fp = fopen(filepath, "rb");
    if (fp == NULL) {
        return NULL;
    }

    if (fseek(fp, 0, SEEK_END) != 0) {
        fclose(fp);
        return NULL;
    }
    long nbytes = ftell(fp);
    if (nbytes < 0) {
        fclose(fp);
        return NULL;
    }
    if (fseek(fp, 0, SEEK_SET) != 0) {
        fclose(fp);
        return NULL;
    }

    char *buf = (char *)malloc((size_t)nbytes + 1);
    if (buf == NULL) {
        fclose(fp);
        return NULL;
    }

    size_t nr = fread(buf, 1, (size_t)nbytes, fp);
    fclose(fp);
    if (nr != (size_t)nbytes) {
        free(buf);
        return NULL;
    }
    buf[nbytes] = '\0';
    return buf;
}

static double get_number_or_default(const cJSON *obj, const char *name, double defval) {
    const cJSON *item = cJSON_GetObjectItemCaseSensitive((cJSON *)obj, name);
    if (cJSON_IsNumber(item)) {
        return item->valuedouble;
    }
    return defval;
}

SimParams parse_params_json(const char *filepath) {
    SimParams p;
    memset(&p, 0, sizeof(p));
    p.tau = INFINITY;

    char *text = read_file_text(filepath);
    if (text == NULL) {
        return p;
    }

    cJSON *root = cJSON_Parse(text);
    free(text);
    if (root == NULL) {
        return p;
    }

    cJSON *params = cJSON_GetObjectItemCaseSensitive(root, "parameters");
    cJSON *sweep = cJSON_GetObjectItemCaseSensitive(root, "bias_sweep");

    p.N = (int)get_number_or_default(params, "N", 0.0);
    p.vmode = get_number_or_default(params, "vmode", 0.0);
    p.alphaL = get_number_or_default(params, "alphaL", 0.0);
    p.alphaR = get_number_or_default(params, "alphaR", 0.0);
    p.lambda = get_number_or_default(params, "lambda", 0.0);
    p.T = get_number_or_default(params, "T", 0.0);
    p.eta = get_number_or_default(params, "eta", 0.0);
    p.Vg = get_number_or_default(params, "Vg", 0.0);

    cJSON *tau = cJSON_GetObjectItemCaseSensitive(params, "tau");
    if (cJSON_IsString(tau) && tau->valuestring != NULL && strcmp(tau->valuestring, "Inf") == 0) {
        p.tau = INFINITY;
    } else if (cJSON_IsNumber(tau)) {
        p.tau = tau->valuedouble;
    }

    p.Vsd_start = get_number_or_default(sweep, "Vsd_start", 0.0);
    p.Vsd_end = get_number_or_default(sweep, "Vsd_end", 0.0);
    p.Vsd_step = get_number_or_default(sweep, "Vsd_step", 0.0);

    cJSON_Delete(root);
    return p;
}

double *load_matlab_vsd(const char *filepath, int *nVsd) {
    if (nVsd != NULL) {
        *nVsd = 0;
    }

    char *text = read_file_text(filepath);
    if (text == NULL) {
        return NULL;
    }

    cJSON *root = cJSON_Parse(text);
    free(text);
    if (root == NULL) {
        return NULL;
    }

    cJSON *arr = cJSON_GetObjectItemCaseSensitive(root, "Vsd");
    if (!cJSON_IsArray(arr)) {
        cJSON_Delete(root);
        return NULL;
    }

    const int n = cJSON_GetArraySize(arr);
    if (n <= 0) {
        cJSON_Delete(root);
        return NULL;
    }

    double *vsd = (double *)malloc((size_t)n * sizeof(double));
    if (vsd == NULL) {
        cJSON_Delete(root);
        return NULL;
    }

    for (int i = 0; i < n; ++i) {
        cJSON *item = cJSON_GetArrayItem(arr, i);
        vsd[i] = cJSON_IsNumber(item) ? item->valuedouble : 0.0;
    }

    cJSON_Delete(root);
    if (nVsd != NULL) {
        *nVsd = n;
    }
    return vsd;
}

void load_matlab_reference(const char *filepath,
                           double **ref_I_tol, double **ref_I_seq,
                           double **ref_I_cot, int *npts) {
    if (ref_I_tol != NULL) {
        *ref_I_tol = NULL;
    }
    if (ref_I_seq != NULL) {
        *ref_I_seq = NULL;
    }
    if (ref_I_cot != NULL) {
        *ref_I_cot = NULL;
    }
    if (npts != NULL) {
        *npts = 0;
    }

    char *text = read_file_text(filepath);
    if (text == NULL) {
        return;
    }

    cJSON *root = cJSON_Parse(text);
    free(text);
    if (root == NULL) {
        return;
    }

    cJSON *a_tol = cJSON_GetObjectItemCaseSensitive(root, "I_tol");
    cJSON *a_seq = cJSON_GetObjectItemCaseSensitive(root, "I_seq");
    cJSON *a_cot = cJSON_GetObjectItemCaseSensitive(root, "I_cot");
    if (!cJSON_IsArray(a_tol) || !cJSON_IsArray(a_seq) || !cJSON_IsArray(a_cot)) {
        cJSON_Delete(root);
        return;
    }

    int n = cJSON_GetArraySize(a_tol);
    int n_seq = cJSON_GetArraySize(a_seq);
    int n_cot = cJSON_GetArraySize(a_cot);
    if (n_seq < n) {
        n = n_seq;
    }
    if (n_cot < n) {
        n = n_cot;
    }
    if (n <= 0) {
        cJSON_Delete(root);
        return;
    }

    double *tol = (double *)malloc((size_t)n * sizeof(double));
    double *seq = (double *)malloc((size_t)n * sizeof(double));
    double *cot = (double *)malloc((size_t)n * sizeof(double));
    if (tol == NULL || seq == NULL || cot == NULL) {
        free(tol);
        free(seq);
        free(cot);
        cJSON_Delete(root);
        return;
    }

    for (int i = 0; i < n; ++i) {
        cJSON *it = cJSON_GetArrayItem(a_tol, i);
        cJSON *is = cJSON_GetArrayItem(a_seq, i);
        cJSON *ic = cJSON_GetArrayItem(a_cot, i);
        tol[i] = cJSON_IsNumber(it) ? it->valuedouble : 0.0;
        seq[i] = cJSON_IsNumber(is) ? is->valuedouble : 0.0;
        cot[i] = cJSON_IsNumber(ic) ? ic->valuedouble : 0.0;
    }

    cJSON_Delete(root);
    if (ref_I_tol != NULL) {
        *ref_I_tol = tol;
    }
    if (ref_I_seq != NULL) {
        *ref_I_seq = seq;
    }
    if (ref_I_cot != NULL) {
        *ref_I_cot = cot;
    }
    if (npts != NULL) {
        *npts = n;
    }
}

static void write_json_array(FILE *fp, const char *name, const double *arr, int n, int comma) {
    fprintf(fp, "  \"%s\": [\n", name);
    for (int i = 0; i < n; ++i) {
        fprintf(fp, "    %.17g%s\n", arr[i], (i + 1 < n) ? "," : "");
    }
    fprintf(fp, "  ]%s\n", comma ? "," : "");
}

void write_results_json(const char *filepath, const char *spec,
                        const SimParams *params, const SimulationResult *res,
                        double wall_time) {
    if (filepath == NULL || spec == NULL || params == NULL || res == NULL) {
        return;
    }

    FILE *fp = fopen(filepath, "w");
    if (fp == NULL) {
        return;
    }

    time_t now = time(NULL);
    struct tm *lt = localtime(&now);
    char tbuf[64] = {0};
    if (lt != NULL) {
        strftime(tbuf, sizeof(tbuf), "%Y-%m-%d %H:%M:%S", lt);
    }

    fprintf(fp, "{\n");
    fprintf(fp, "  \"spec\": \"%s\",\n", spec);
    fprintf(fp, "  \"language\": \"C\",\n");
    fprintf(fp, "  \"version\": \"gcc %s\",\n", __VERSION__);
    fprintf(fp, "  \"wall_time_seconds\": %.17g,\n", wall_time);
    fprintf(fp, "  \"parameters\": {\n");
    fprintf(fp, "    \"N\": %d,\n", params->N);
    fprintf(fp, "    \"vmode\": %.17g,\n", params->vmode);
    fprintf(fp, "    \"alphaL\": %.17g,\n", params->alphaL);
    fprintf(fp, "    \"alphaR\": %.17g,\n", params->alphaR);
    fprintf(fp, "    \"lambda\": %.17g,\n", params->lambda);
    fprintf(fp, "    \"T\": %.17g,\n", params->T);
    fprintf(fp, "    \"eta\": %.17g,\n", params->eta);
    fprintf(fp, "    \"Vg\": %.17g,\n", params->Vg);
    if (isinf(params->tau)) {
        fprintf(fp, "    \"tau\": \"Inf\"\n");
    } else {
        fprintf(fp, "    \"tau\": %.17g\n", params->tau);
    }
    fprintf(fp, "  },\n");
    fprintf(fp, "  \"bias_sweep\": {\n");
    fprintf(fp, "    \"Vsd_start\": %.17g,\n", params->Vsd_start);
    fprintf(fp, "    \"Vsd_end\": %.17g,\n", params->Vsd_end);
    fprintf(fp, "    \"Vsd_step\": %.17g\n", params->Vsd_step);
    fprintf(fp, "  },\n");

    write_json_array(fp, "Vsd", res->Vsd, res->nVsd, 1);
    write_json_array(fp, "I_tol", res->I_tol, res->nVsd, 1);
    write_json_array(fp, "I_seq", res->I_seq, res->nVsd, 1);
    write_json_array(fp, "I_cot", res->I_cot, res->nVsd, 1);
    fprintf(fp, "  \"timestamp\": \"%s\"\n", tbuf);
    fprintf(fp, "}\n");

    fclose(fp);
}

void write_results_csv(const char *filepath, const SimulationResult *res) {
    if (filepath == NULL || res == NULL) {
        return;
    }

    FILE *fp = fopen(filepath, "w");
    if (fp == NULL) {
        return;
    }

    fprintf(fp, "Vsd_V,I_tol_A,I_seq_A,I_cot_A\n");
    for (int i = 0; i < res->nVsd; ++i) {
        fprintf(fp, "%.6e,%.6e,%.6e,%.6e\n",
                res->Vsd[i], res->I_tol[i], res->I_seq[i], res->I_cot[i]);
    }

    fclose(fp);
}
