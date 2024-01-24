#!/bin/bash
# Compara caracteriticas especificas entre estos archivos (ver repositorio)
# Definir la ruta de los archivos
ruta_archivos="/Escritorio/structures/"

# Definir el archivo de referencia
archivo_referencia="whitout_ter.pdb"

# Definir los archivos a comparar
archivos_comparar=("I1_with_opsPEC.pdb" "I2_with_opsPEC.pdb" "I3_with_opsPEC.pdb" "I4_with_opsPEC.pdb")

# Definir las características a unificar
caracteristicas="ATOM"

# Iterar sobre los archivos a comparar y unificar las características
for archivo in "${archivos_comparar[@]}"
do
    # Crear un nuevo archivo con las características unificadas
    cat "${ruta_archivos}${archivo}" | grep "${caracteristicas}" > "${ruta_archivos}${archivo}_unificado.pdb"

    # Comparar el archivo unificado con el archivo de referencia utilizando el programa pdbalign
    pdbalign -s "${ruta_archivos}${archivo}_unificado.pdb" -r "${ruta_archivos}${archivo_referencia}" -o "${ruta_archivos}${archivo}_alineado.pdb"
done
