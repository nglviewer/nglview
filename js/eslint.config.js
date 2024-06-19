module.exports = {
  parser: '@typescript-eslint/parser',
  extends: 'eslint:recommended',
  rules: {
    'no-console': 'off',
    'no-constant-condition': 'off',
    'no-unused-vars': 'off',
    'no-redeclare': 'off',
    'valid-jsdoc': 'error',
  },
  env: {
    es6: true,
    node: true,
    mocha: true,
    browser: true,
    worker: true,
  },
  parserOptions: {
    sourceType: 'module',
  },
  ignorePatterns: ['node_modules', 'dist', 'coverage', '**/*.d.ts'],
};
