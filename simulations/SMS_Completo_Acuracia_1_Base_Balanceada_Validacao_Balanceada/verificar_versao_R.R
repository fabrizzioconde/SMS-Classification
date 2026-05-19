# Script para verificar a versão do R
cat("=== Informações da Versão do R ===\n\n")
cat("Versão completa:\n")
print(R.version.string)
cat("\nVersão detalhada:\n")
print(R.version)
cat("\nVersão major.minor:\n")
cat(paste(R.version$major, R.version$minor, sep = "."), "\n")
