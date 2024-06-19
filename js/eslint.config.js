module.exports = {
  files: [
      "**/*.ts",
      "**/*.cts",
      "**.*.mts"
  ],
  rules: {
    'no-console': 'off',
    'no-constant-condition': 'off',
    'no-unused-vars': 'off',
    'no-redeclare': 'off',
    'valid-jsdoc': 'error',
  },
  ignorePatterns: ['node_modules', 'dist', 'coverage']

};
