import numpy as np
import matplotlib.pyplot as plt  # por compatibilidad, aunque usamos Figure
import pandas as pd
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure

# -----------------------------
# Constantes y nombres de meses
# -----------------------------
MESES = [
    "Enero", "Febrero", "Marzo", "Abril",
    "Mayo", "Junio", "Julio", "Agosto",
    "Septiembre", "Octubre", "Noviembre", "Diciembre"
]
DAYS_POR_MES = 30  # aproximación para mostrar meses en el eje X


# ===========================
# FUNCIONES DEL MODELO
# ===========================

def beta_sazonal(t, beta0, fuerza_estacional):
    """
    Función de transmisión estacional.
    t: tiempo en días
    beta0: tasa base de transmisión
    fuerza_estacional: amplitud (0 a 1)
    """
    return beta0 * (1.0 + fuerza_estacional * np.sin(2.0 * np.pi * t / 365.0))


def correr_simulacion(
    N, I0, E0, R0,
    beta0,
    dias_incubacion, dias_infecciosos,
    dias,
    fuerza_estacional,
    anio=None,
    nombre=None,
    days_por_mes=DAYS_POR_MES,
):
    """
    Modelo SEIR con integración por Euler (paso 1 día).
    NO dibuja gráficas, sólo devuelve los resultados y estadísticas.
    """
    try:
        # Validaciones
        if N <= 0:
            messagebox.showwarning("Dato inválido", "La población total debe ser mayor que 0.")
            return None
        if dias <= 0:
            messagebox.showwarning("Dato inválido", "Los días de simulación deben ser mayores que 0.")
            return None
        if dias_incubacion <= 0 or dias_infecciosos <= 0:
            messagebox.showwarning(
                "Dato inválido",
                "Los días de incubación y de etapa contagiosa deben ser mayores que 0."
            )
            return None

        sigma = 1.0 / dias_incubacion
        gamma = 1.0 / dias_infecciosos
        dias = int(dias)

        S0 = N - I0 - E0 - R0
        if S0 < 0:
            messagebox.showwarning(
                "Dato inválido",
                "La suma de infectados, expuestos e inmunes iniciales es mayor que la población total."
            )
            return None

        # Arrays
        t = np.arange(0, dias)
        S = np.zeros(dias)
        E = np.zeros(dias)
        I = np.zeros(dias)
        R = np.zeros(dias)
        beta_t = np.zeros(dias)

        S[0] = S0
        E[0] = E0
        I[0] = I0
        R[0] = R0
        beta_t[0] = beta_sazonal(0, beta0, fuerza_estacional)

        # Integración por Euler
        for day in range(1, dias):
            beta_val = beta_sazonal(day, beta0, fuerza_estacional)
            beta_t[day] = beta_val

            dS = -beta_val * S[day - 1] * I[day - 1] / N
            dE = beta_val * S[day - 1] * I[day - 1] / N - sigma * E[day - 1]
            dI = sigma * E[day - 1] - gamma * I[day - 1]
            dR = gamma * I[day - 1]

            S[day] = S[day - 1] + dS
            E[day] = E[day - 1] + dE
            I[day] = I[day - 1] + dI
            R[day] = R[day - 1] + dR

        # Estadísticas de pico
        if np.all(I == 0):
            peak_idx = 0
        else:
            peak_idx = int(np.argmax(I))  # índice del día del pico

        peak_day = peak_idx
        peak_value = float(I[peak_idx])
        peak_month_idx = int((peak_day // days_por_mes) % 12)  # 0-11
        peak_month_name = MESES[peak_month_idx]

        total_recuperados = float(R[-1])
        final_infectados = float(I[-1])
        total_casos_estimados = total_recuperados + final_infectados

        # --------------------------
        # Construcción de títulos
        # --------------------------
        nombre_str = nombre if (nombre is not None and str(nombre).strip() != "") else "escenario"

        rango_str = ""
        approx_years = None

        if anio is not None:
            anio_base = int(anio)         # año de los datos (ej. 2024)
            anio_inicio = anio_base + 1   # año en que empieza la simulación (ej. 2025)

            # Años aproximados que cubre la simulación
            duracion_anios = dias / 365.0
            approx_years = max(1, int(np.ceil(duracion_anios)))
            anio_fin = anio_inicio + approx_years - 1

            if approx_years == 1:
                # solo un año (ej. 365 días o menos)
                rango_str = f"{anio_inicio}"
            else:
                # rango aproximado (ej. 1000 días -> ~3 años)
                rango_str = f"{anio_inicio}-{anio_fin}"

        # Elegimos el texto del título
        if rango_str:
            # hay información de años
            es_anio_completo = (dias >= 360 and dias <= 370 and approx_years == 1)
            if es_anio_completo:
                titulo_seir = f"Modelo estacional {nombre_str} {rango_str}"
                titulo_beta = f"Tasa de transmisión β(t) - {nombre_str} {rango_str}"
            else:
                titulo_seir = f"Modelo SEIR {nombre_str} {rango_str} (~{dias} días)"
                titulo_beta = f"Tasa de transmisión β(t) - {nombre_str} {rango_str} (~{dias} días)"
        else:
            # sin año de referencia
            if dias >= 360 and dias <= 370:
                titulo_seir = "Modelo estacional SEIR"
                titulo_beta = "Tasa de transmisión β(t) estacional"
            else:
                titulo_seir = f"Modelo SEIR (~{dias} días)"
                titulo_beta = f"Tasa de transmisión β(t) (~{dias} días)"

        resultados = {
            "t": t,
            "S": S,
            "E": E,
            "I": I,
            "R": R,
            "beta_t": beta_t,
            "peak_day": peak_day,
            "peak_value": peak_value,
            "peak_month_name": peak_month_name,
            "total_recuperados": total_recuperados,
            "final_infectados": final_infectados,
            "total_casos_estimados": total_casos_estimados,
            "dias": dias,
            "anio": anio,
            "nombre": nombre,
            "titulo_seir": titulo_seir,
            "titulo_beta": titulo_beta,
        }
        return resultados

    except Exception as e:
        messagebox.showerror("Error en simulación", f"Ocurrió un error:\n{e}")
        return None


# ===========================
# FUNCIONES DE LA GUI
# ===========================

last_results = None   # último resultado de simulación
current_plot_mode = None  # "seir" o "beta"
fig = None
ax = None
canvas = None
hover_cursor = None   # cursor de mplcursors para hover


def leer_float(entry, nombre_campo):
    """
    Intenta leer un float desde una Entry.
    Lanza ValueError si está vacío o no es número.
    """
    texto = entry.get().strip()
    if texto == "":
        raise ValueError(f"El campo '{nombre_campo}' está vacío.")
    try:
        return float(texto)
    except ValueError:
        raise ValueError(f"El campo '{nombre_campo}' debe ser numérico.")


def leer_int_opcional(entry):
    """
    Devuelve un int o None si la caja está vacía.
    """
    texto = entry.get().strip()
    if texto == "":
        return None
    try:
        return int(float(texto))
    except ValueError:
        messagebox.showerror("Error de datos", "El campo de año debe ser numérico.")
        return None


def actualizar_labels_resultados(res):
    label_peak_val.config(text=f"Pico de personas enfermas: {res['peak_value']:.0f}")
    label_peak_mes.config(text=f"Mes del pico: {res['peak_month_name']}")
    label_final_rec.config(text=f"Personas recuperadas al final: {res['total_recuperados']:.0f}")
    label_total_casos.config(text=f"Casos acumulados aproximados: {res['total_casos_estimados']:.0f}")


def preparar_eje_meses(ax_local, dias):
    """
    Configura el eje X en función de meses con nombres reales.
    """
    meses_totales = int(np.ceil(dias / float(DAYS_POR_MES)))
    max_ticks = 24
    step_meses = 1
    if meses_totales > max_ticks:
        step_meses = int(np.ceil(meses_totales / max_ticks))

    month_positions = np.arange(0, meses_totales, step_meses) * DAYS_POR_MES
    month_positions = month_positions[month_positions < dias]

    month_indices = ((month_positions // DAYS_POR_MES) % 12).astype(int)
    month_labels = [MESES[i] for i in month_indices]

    ax_local.set_xticks(month_positions)
    ax_local.set_xticklabels(month_labels, rotation=45, ha="right")
    ax_local.set_xlabel("Mes")


def limpiar_hover():
    """
    Elimina el cursor/hover actual (si existe) y sus anotaciones.
    """
    global hover_cursor
    if hover_cursor is not None:
        try:
            hover_cursor.remove()
        except Exception:
            pass
        hover_cursor = None


def plot_seir():
    global last_results, ax, canvas, current_plot_mode, hover_cursor

    if last_results is None:
        messagebox.showinfo("Sin datos", "Primero ejecuta una simulación.")
        return

    # limpiar cualquier hover anterior
    limpiar_hover()

    current_plot_mode = "seir"
    res = last_results
    t = res["t"]
    S = res["S"]
    E = res["E"]
    I = res["I"]
    R = res["R"]
    dias = res["dias"]

    ax.clear()

    lnS, = ax.plot(t, S, label="Personas susceptibles (S)")
    lnE, = ax.plot(t, E, label="Personas en incubación (E)")
    lnI, = ax.plot(t, I, label="Personas enfermas (I)")
    lnR, = ax.plot(t, R, label="Personas inmunes/recuperadas (R)")

    preparar_eje_meses(ax, dias)
    ax.set_ylabel("Número de personas")
    ax.set_title(res["titulo_seir"])
    ax.legend()
    ax.grid(True, alpha=0.3)

    # marcar pico
    peak_day = res["peak_day"]
    peak_value = res["peak_value"]
    peak_month_name = res["peak_month_name"]

    ax.plot(peak_day, peak_value, marker='o', markersize=6, color='red')
    ax.axvline(peak_day, linestyle='--', alpha=0.6)
    text_x = min(peak_day + DAYS_POR_MES * 0.2, dias - 1)
    ax.annotate(
        f"Pico: {peak_value:.0f}\n{peak_month_name}",
        xy=(peak_day, peak_value),
        xytext=(text_x, peak_value),
        arrowprops=dict(arrowstyle="->", alpha=0.7),
        bbox=dict(boxstyle="round,pad=0.3", alpha=0.2)
    )

    canvas.draw()

    # Interactividad opcional
    if var_interactive.get():
        try:
            import mplcursors
            hover_cursor = mplcursors.cursor([lnS, lnE, lnI, lnR], hover=True)
            hover_cursor.connect(
                "add",
                lambda sel: sel.annotation.set_text(
                    f"{sel.artist.get_label()}\nDía {int(sel.target[0])}\nValor {sel.target[1]:.0f}"
                )
            )
        except Exception:
            messagebox.showinfo(
                "Interactividad no disponible",
                "Para activar hover instala la librería 'mplcursors' (pip install mplcursors)."
            )


def plot_beta():
    global last_results, ax, canvas, current_plot_mode, hover_cursor

    if last_results is None:
        messagebox.showinfo("Sin datos", "Primero ejecuta una simulación.")
        return

    # limpiar cualquier hover anterior
    limpiar_hover()

    current_plot_mode = "beta"
    res = last_results
    t = res["t"]
    beta_t = res["beta_t"]
    dias = res["dias"]

    ax.clear()

    lnB, = ax.plot(t, beta_t, label="Tasa estacional de transmisión")
    preparar_eje_meses(ax, dias)
    ax.set_ylabel("Tasa de transmisión β(t) (1/día)")
    ax.set_title(res["titulo_beta"])
    ax.grid(True, alpha=0.3)
    ax.legend()

    canvas.draw()

    if var_interactive.get():
        try:
            import mplcursors
            hover_cursor = mplcursors.cursor([lnB], hover=True)
            hover_cursor.connect(
                "add",
                lambda sel: sel.annotation.set_text(
                    f"Tasa de transmisión\nDía {int(sel.target[0])}\nValor {sel.target[1]:.3f} 1/día"
                )
            )
        except Exception:
            messagebox.showinfo(
                "Interactividad no disponible",
                "Para activar hover instala la librería 'mplcursors' (pip install mplcursors)."
            )


def on_hover_toggle(*args):
    """
    Cuando el usuario marca/desmarca el checkbox de interactividad,
    redibuja la última gráfica respetando el nuevo estado.
    """
    if last_results is None or current_plot_mode is None:
        return

    if current_plot_mode == "seir":
        plot_seir()
    elif current_plot_mode == "beta":
        plot_beta()


def simular_desde_campos():
    """
    Lee todos los campos de texto y corre la simulación.
    """
    global last_results

    try:
        nombre = entry_nombre.get().strip()
        anio = leer_int_opcional(entry_anio)
        if anio is None and entry_anio.get().strip() != "":
            # ya se mostró el error, abortamos
            return

        N = leer_float(entry_N, "Población total")
        I0 = leer_float(entry_I0, "Personas enfermas al inicio")
        E0 = leer_float(entry_E0, "Personas en incubación al inicio")
        R0 = leer_float(entry_R0, "Personas inmunes al inicio")
        beta0 = leer_float(entry_beta0, "Intensidad de contagio base (β0)")
        dias_incubacion = leer_float(entry_dias_incubacion, "Días de incubación")
        dias_infecciosos = leer_float(entry_dias_infecciosos, "Días de etapa contagiosa")
        dias = leer_float(entry_dias, "Duración de la simulación (días)")
        fuerza_estacional = leer_float(entry_fuerza, "Fuerza del efecto climático")

        res = correr_simulacion(
            N, I0, E0, R0,
            beta0,
            dias_incubacion, dias_infecciosos,
            dias,
            fuerza_estacional,
            anio=anio,
            nombre=nombre if nombre != "" else None,
            days_por_mes=DAYS_POR_MES,
        )

        if res:
            last_results = res
            actualizar_labels_resultados(res)
            plot_seir()  # mostramos por defecto S/E/I/R

    except ValueError as e:
        messagebox.showerror("Error de datos", str(e))


def cargar_excel_y_simular():
    """
    Lee parámetros desde un archivo Excel o CSV y corre la simulación
    usando la PRIMERA FILA.
    NO rellena las cajas de texto: los campos quedan para edición manual.
    - Si es Excel: lee la hoja llamada 'Datos'.
    - Si es CSV: lee el archivo directo.
    """
    global last_results

    ruta = filedialog.askopenfilename(
        title="Selecciona archivo de Excel/CSV con los parámetros",
        filetypes=[("Archivos de Excel/CSV", "*.xlsx *.xls *.csv")]
    )
    if not ruta:
        return  # usuario canceló

    try:
        if ruta.lower().endswith(".csv"):
            df = pd.read_csv(ruta)
        else:
            # Leemos hoja "Datos"
            try:
                df = pd.read_excel(ruta, sheet_name="Datos")
            except ValueError:
                messagebox.showerror(
                    "Error en la hoja",
                    "El archivo de Excel debe tener una hoja llamada 'Datos'."
                )
                return

        columnas_esperadas = [
            "Año del escenario", "Nombre del lugar",
            "Población total (N)", "Infectados iniciales (I0)",
            "Expuestos iniciales (E0)", "Recuperados iniciales (R0)",
            "Tasa base de transmisión (beta0)",
            "Días de incubación", "Días infecciosos",
            "Días de simulación",
            "Fuerza estacional (0 a 1)",
        ]
        for col in columnas_esperadas:
            if col not in df.columns:
                messagebox.showerror(
                    "Columnas faltantes",
                    "La hoja 'Datos' debe tener estas columnas:\n" +
                    ", ".join(columnas_esperadas)
                )
                return

        if len(df) == 0:
            messagebox.showerror("Archivo vacío", "La hoja 'Datos' no tiene filas de datos.")
            return

        fila = df.iloc[0]

        anio = int(fila["Año del escenario"])
        nombre = str(fila["Nombre del lugar"])

        N = float(fila["Población total (N)"])
        I0 = float(fila["Infectados iniciales (I0)"])
        E0 = float(fila["Expuestos iniciales (E0)"])
        R0 = float(fila["Recuperados iniciales (R0)"])
        beta0 = float(fila["Tasa base de transmisión (beta0)"])
        dias_incubacion = float(fila["Días de incubación"])
        dias_infecciosos = float(fila["Días infecciosos"])
        dias = float(fila["Días de simulación"])
        fuerza_estacional = float(fila["Fuerza estacional (0 a 1)"])

        # Ejecutar simulación directamente (sin rellenar textfields)
        res = correr_simulacion(
            N, I0, E0, R0,
            beta0,
            dias_incubacion, dias_infecciosos,
            dias,
            fuerza_estacional,
            anio=anio,
            nombre=nombre,
            days_por_mes=DAYS_POR_MES,
        )

        if res:
            last_results = res
            actualizar_labels_resultados(res)
            plot_seir()

    except Exception as e:
        messagebox.showerror("Error al leer archivo", f"No se pudo leer el archivo:\n{e}")


def generar_formato_vacio():
    """
    Genera un archivo Excel con DOS hojas:
    - 'Instrucciones': explicación de cada campo.
    - 'Datos': plantilla vacía con las columnas a rellenar (nombres descriptivos).
    """
    columnas_datos = [
        "Año del escenario", "Nombre del lugar",
        "Población total (N)", "Infectados iniciales (I0)",
        "Expuestos iniciales (E0)", "Recuperados iniciales (R0)",
        "Tasa base de transmisión (beta0)",
        "Días de incubación", "Días infecciosos",
        "Días de simulación",
        "Fuerza estacional (0 a 1)",
    ]

    # Hoja de instrucciones
    filas_instr = [
        {
            "Campo": "Año del escenario",
            "Descripción": "Año de referencia del escenario que deseas simular (por ejemplo 2024)."
        },
        {
            "Campo": "Nombre del lugar",
            "Descripción": "Nombre del lugar o región (por ejemplo 'Brasil', 'Oaxaca de Juárez')."
        },
        {
            "Campo": "Población total (N)",
            "Descripción": "Número total de habitantes del lugar."
        },
        {
            "Campo": "Infectados iniciales (I0)",
            "Descripción": "Personas enfermas al inicio de la simulación (día 0)."
        },
        {
            "Campo": "Expuestos iniciales (E0)",
            "Descripción": "Personas en incubación (infectadas pero aún sin síntomas) al inicio."
        },
        {
            "Campo": "Recuperados iniciales (R0)",
            "Descripción": "Personas que ya son inmunes o se han recuperado antes del inicio."
        },
        {
            "Campo": "Tasa base de transmisión (beta0)",
            "Descripción": "Intensidad promedio del contagio por día (β0)."
        },
        {
            "Campo": "Días de incubación",
            "Descripción": "Duración promedio de la incubación de la enfermedad (en días)."
        },
        {
            "Campo": "Días infecciosos",
            "Descripción": "Duración promedio de la etapa en la que la persona contagia (en días)."
        },
        {
            "Campo": "Días de simulación",
            "Descripción": "Cuántos días totales quieres simular (ej. 365 para un año)."
        },
        {
            "Campo": "Fuerza estacional (0 a 1)",
            "Descripción": "Qué tan fuerte afecta el clima a los contagios (0 = nada, 1 = muy fuerte)."
        },
        {
            "Campo": "",
            "Descripción": (
                "IMPORTANTE: la simulación leerá SIEMPRE la hoja 'Datos' "
                "y usará solamente la PRIMERA fila de parámetros que encuentre."
            ),
        },
    ]
    df_instr = pd.DataFrame(filas_instr, columns=["Campo", "Descripción"])

    # Hoja de datos vacía
    df_datos = pd.DataFrame(columns=columnas_datos)

    ruta = filedialog.asksaveasfilename(
        title="Guardar plantilla de parámetros",
        defaultextension=".xlsx",
        filetypes=[("Archivo de Excel", "*.xlsx")]
    )
    if not ruta:
        return  # usuario canceló

    try:
        with pd.ExcelWriter(ruta, engine="openpyxl") as writer:
            df_instr.to_excel(writer, sheet_name="Instrucciones", index=False)
            df_datos.to_excel(writer, sheet_name="Datos", index=False)

        messagebox.showinfo(
            "Plantilla creada",
            f"Se creó una plantilla con las hojas 'Instrucciones' y 'Datos':\n{ruta}"
        )
    except Exception as e:
        messagebox.showerror("Error al guardar plantilla", f"No se pudo guardar el archivo:\n{e}")


# ===========================
# GUI: VENTANA PRINCIPAL
# ===========================

root = tk.Tk()
root.title("Simulador SEIR de Dengue - Estacional")

# Maximizar ventana al inicio
try:
    root.state("zoomed")  # Windows
except tk.TclError:
    try:
        root.attributes("-zoomed", True)  # Linux
    except tk.TclError:
        pass  # si no se puede, la dejamos normal

# --------- Frame principal del simulador (aún sin mostrar) ---------
main_frame = ttk.Frame(root, padding=10)

# Dividimos en dos columnas: izquierda parámetros, derecha gráficas
left_frame = ttk.Frame(main_frame)
left_frame.grid(row=0, column=0, sticky="nsw", padx=(0, 10))

right_frame = ttk.Frame(main_frame)
right_frame.grid(row=0, column=1, sticky="nsew")

main_frame.columnconfigure(1, weight=1)
main_frame.rowconfigure(0, weight=1)

# --------- LADO IZQUIERDO: PARÁMETROS Y BOTONES ---------
row = 0
ttk.Label(left_frame, text="Simulador SEIR de Dengue", font=("Arial", 12, "bold")).grid(
    row=row, column=0, columnspan=2, pady=(0, 10)
)

row += 1
ttk.Label(left_frame, text="Lugar o escenario (ej. Brasil):").grid(row=row, column=0, sticky="e", pady=2)
entry_nombre = ttk.Entry(left_frame)
entry_nombre.grid(row=row, column=1, sticky="w", pady=2)

row += 1
ttk.Label(left_frame, text="Año del escenario (ej. 2024):").grid(row=row, column=0, sticky="e", pady=2)
entry_anio = ttk.Entry(left_frame)
entry_anio.grid(row=row, column=1, sticky="w", pady=2)

row += 1
ttk.Label(left_frame, text="Población total (personas):").grid(row=row, column=0, sticky="e", pady=2)
entry_N = ttk.Entry(left_frame)
entry_N.grid(row=row, column=1, sticky="w", pady=2)

row += 1
ttk.Label(left_frame, text="Personas enfermas al inicio:").grid(row=row, column=0, sticky="e", pady=2)
entry_I0 = ttk.Entry(left_frame)
entry_I0.grid(row=row, column=1, sticky="w", pady=2)

row += 1
ttk.Label(left_frame, text="Personas en incubación al inicio:").grid(row=row, column=0, sticky="e", pady=2)
entry_E0 = ttk.Entry(left_frame)
entry_E0.grid(row=row, column=1, sticky="w", pady=2)

row += 1
ttk.Label(left_frame, text="Personas inmunes al inicio:").grid(row=row, column=0, sticky="e", pady=2)
entry_R0 = ttk.Entry(left_frame)
entry_R0.grid(row=row, column=1, sticky="w", pady=2)

row += 1
ttk.Label(left_frame, text="Intensidad de contagio base (β0):").grid(row=row, column=0, sticky="e", pady=2)
entry_beta0 = ttk.Entry(left_frame)
entry_beta0.grid(row=row, column=1, sticky="w", pady=2)

row += 1
ttk.Label(left_frame, text="Días de incubación:").grid(row=row, column=0, sticky="e", pady=2)
entry_dias_incubacion = ttk.Entry(left_frame)
entry_dias_incubacion.grid(row=row, column=1, sticky="w", pady=2)

row += 1
ttk.Label(left_frame, text="Días de etapa contagiosa:").grid(row=row, column=0, sticky="e", pady=2)
entry_dias_infecciosos = ttk.Entry(left_frame)
entry_dias_infecciosos.grid(row=row, column=1, sticky="w", pady=2)

row += 1
ttk.Label(left_frame, text="Duración de la simulación (días):").grid(row=row, column=0, sticky="e", pady=2)
entry_dias = ttk.Entry(left_frame)
entry_dias.grid(row=row, column=1, sticky="w", pady=2)

row += 1
ttk.Label(left_frame, text="Fuerza del efecto climático (0 a 1):").grid(row=row, column=0, sticky="e", pady=2)
entry_fuerza = ttk.Entry(left_frame)
entry_fuerza.grid(row=row, column=1, sticky="w", pady=2)

row += 1
btn_simular = ttk.Button(left_frame, text="Simular usando los datos de la izquierda",
                         command=simular_desde_campos)
btn_simular.grid(row=row, column=0, columnspan=2, pady=(10, 5), sticky="ew")

row += 1
btn_excel = ttk.Button(left_frame, text="Cargar archivo Excel/CSV y simular",
                       command=cargar_excel_y_simular)
btn_excel.grid(row=row, column=0, columnspan=2, pady=5, sticky="ew")

row += 1
btn_formato = ttk.Button(left_frame, text="Crear plantilla de Excel para rellenar",
                         command=generar_formato_vacio)
btn_formato.grid(row=row, column=0, columnspan=2, pady=5, sticky="ew")

# Interactividad (hover)
var_interactive = tk.BooleanVar(value=False)
var_interactive.trace_add("write", on_hover_toggle)

row += 1
chk_hover = ttk.Checkbutton(
    left_frame,
    text="Mostrar información al pasar el mouse sobre la gráfica",
    variable=var_interactive
)
chk_hover.grid(row=row, column=0, columnspan=2, pady=(5, 10), sticky="w")

# Etiquetas de resultados
row += 1
ttk.Label(left_frame, text="Resumen de la simulación:", font=("Arial", 11, "bold")).grid(
    row=row, column=0, columnspan=2, sticky="w"
)

row += 1
label_peak_val = ttk.Label(left_frame, text="Pico de personas enfermas: —")
label_peak_val.grid(row=row, column=0, sticky="w", pady=2)
label_peak_mes = ttk.Label(left_frame, text="Mes del pico: —")
label_peak_mes.grid(row=row, column=1, sticky="w", pady=2)

row += 1
label_final_rec = ttk.Label(left_frame, text="Personas recuperadas al final: —")
label_final_rec.grid(row=row, column=0, sticky="w", pady=2)
label_total_casos = ttk.Label(left_frame, text="Casos acumulados aproximados: —")
label_total_casos.grid(row=row, column=1, sticky="w", pady=2)

left_frame.columnconfigure(0, weight=0)
left_frame.columnconfigure(1, weight=1)

# --------- LADO DERECHO: GRÁFICAS ---------
ttk.Label(right_frame, text="Gráficas de la simulación", font=("Arial", 12, "bold")).pack(pady=(0, 5))

fig = Figure(figsize=(6, 5), dpi=100)
ax = fig.add_subplot(111)
ax.text(0.5, 0.5, "Ejecuta una simulación para ver las gráficas",
        ha="center", va="center", transform=ax.transAxes)
ax.set_axis_off()

canvas = FigureCanvasTkAgg(fig, master=right_frame)
canvas_widget = canvas.get_tk_widget()
canvas_widget.pack(fill="both", expand=True)
canvas.draw()

buttons_frame = ttk.Frame(right_frame)
buttons_frame.pack(fill="x", pady=(5, 0))

btn_plot_seir = ttk.Button(
    buttons_frame,
    text="Ver número de personas en cada grupo (S/E/I/R)",
    command=plot_seir
)
btn_plot_seir.pack(side="left", expand=True, fill="x", padx=(0, 5))

btn_plot_beta = ttk.Button(
    buttons_frame,
    text="Ver tasa estacional de transmisión",
    command=plot_beta
)
btn_plot_beta.pack(side="left", expand=True, fill="x", padx=(5, 0))

right_frame.columnconfigure(0, weight=1)


# ===========================
# PANTALLA DE BIENVENIDA
# ===========================

def ir_a_principal():
    welcome_frame.pack_forget()
    main_frame.pack(fill="both", expand=True)


welcome_frame = ttk.Frame(root, padding=20)
welcome_frame.pack(fill="both", expand=True)

texto_bienvenida = (
    "Bienvenido al Simulador SEIR de Dengue\n\n"
    "Este programa permite explorar cómo podría propagarse el dengue en un lugar y año determinados.\n\n"
    "El modelo divide a la población en cuatro grupos:\n"
    "  • S: personas susceptibles (aún sanas)\n"
    "  • E: personas en incubación (infectadas, sin síntomas)\n"
    "  • I: personas enfermas (pueden contagiar)\n"
    "  • R: personas inmunes o recuperadas\n\n"
    "Puedes:\n"
    "  1) Capturar los datos manualmente (población, casos iniciales, duración, efecto del clima), o\n"
    "  2) Cargar un archivo de Excel con esos valores usando la plantilla incluida.\n\n"
    "El simulador genera dos gráficas:\n"
    "  • Número de personas en cada grupo (S/E/I/R) a lo largo del tiempo.\n"
    "  • La tasa de transmisión β(t), que indica qué tan fácil se contagia el dengue según el clima.\n\n"
    "Usa este programa sólo como herramienta educativa para comparar escenarios, "
    "no como predicción exacta."
)

label_bienvenida = ttk.Label(
    welcome_frame,
    text=texto_bienvenida,
    justify="left",
    font=("Arial", 11)
)
label_bienvenida.pack(pady=10)

btn_inicio = ttk.Button(welcome_frame, text="Iniciar simulador", command=ir_a_principal)
btn_inicio.pack(pady=10)

btn_salir = ttk.Button(welcome_frame, text="Salir", command=root.destroy)
btn_salir.pack(pady=5)

root.mainloop()